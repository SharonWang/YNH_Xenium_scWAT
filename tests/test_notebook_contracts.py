import json
import importlib.util
import re
import subprocess
import sys
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
NOTEBOOKS = REPO / "notebooks"


def read_notebook(path):
    return json.loads(path.read_text(encoding="utf-8"))


def parameter_source(notebook):
    cells = [
        cell
        for cell in notebook["cells"]
        if "parameters" in cell.get("metadata", {}).get("tags", [])
    ]
    if len(cells) != 1:
        raise AssertionError("exactly one parameter cell is required")
    return "".join(cells[0]["source"])


def notebook_text(notebook, cell_type=None):
    return "\n".join(
        "".join(cell.get("source", []))
        for cell in notebook["cells"]
        if cell_type is None or cell.get("cell_type") == cell_type
    )


class NotebookContracts(unittest.TestCase):
    def test_renderer_templates_share_the_fixed_qc_contract(self):
        renderer_path = REPO / "scripts" / "render_notebooks.py"
        spec = importlib.util.spec_from_file_location("scwat_renderer", renderer_path)
        renderer = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(renderer)

        section_text = notebook_text(renderer.notebook(renderer.section_cells(False)))
        for token in [
            "FIXED_CELL_QC_THRESHOLDS_PATH",
            "read_fixed_cell_qc_thresholds",
            "fixed_thresholds = fixed_thresholds",
            "5 < nFeature_Xenium < 200",
            "10 < nCount_Xenium < 1000",
        ]:
            self.assertIn(token, section_text)

        summary_text = notebook_text(renderer.notebook(renderer.summary_cells(False)))
        self.assertIn("FIXED_CELL_QC_THRESHOLDS_PATH", summary_text)
        self.assertIn("apply_fixed_primary_bounds", summary_text)
        self.assertNotIn("build_cell_downstream_masks(", summary_text)

    def test_four_region_notebooks_are_parameter_locked(self):
        expected = {f"Region_{index}" for index in range(1, 5)}
        observed = set()

        for index in range(1, 5):
            path = NOTEBOOKS / f"01_QC_Region{index}.ipynb"
            self.assertTrue(path.exists(), f"missing committed notebook: {path.name}")
            notebook = read_notebook(path)
            match = re.search(
                r'^REGION_ID\s*<-\s*"(Region_[1-4])"$',
                parameter_source(notebook),
                re.MULTILINE,
            )
            self.assertIsNotNone(match, f"fixed REGION_ID missing from {path.name}")
            observed.add(match.group(1))
            self.assertTrue(
                notebook["metadata"]["kernelspec"]["name"].startswith("ir"),
                "an R kernelspec is required",
            )
            text = notebook_text(notebook)
            code = notebook_text(notebook, "code")
            for token in [
                "read_fixed_cell_qc_thresholds",
                "fixed_cell_qc_thresholds.tsv",
                "fixed_thresholds = fixed_thresholds",
                "transcript_gene_quality <- data.frame()",
                "build_cell_downstream_masks",
                "cell_downstream_masks.tsv.gz",
                "section_downstream_decision.tsv",
                "raw objects remain unchanged",
                "5 < nFeature_Xenium < 200",
                "10 < nCount_Xenium < 1000",
            ]:
                self.assertIn(token, text)
            self.assertNotRegex(code, r"(?m)^primary_include\s*<-")

        self.assertEqual(observed, expected)

    def test_summary_notebook_is_the_complete_reader_facing_report(self):
        notebook = read_notebook(NOTEBOOKS / "02_slide_QC_summary.ipynb")
        text = notebook_text(notebook)
        code = notebook_text(notebook, "code")
        required_headings = [
            "## Setup",
            "## Step1 - Inputs and validation",
            "## Step2 - QC decision",
            "## Step3 - Alarm evidence and evidence-only gene tiers",
            "## Step 4 - Subset versus full-data burden",
            "## Step5 - Spatial QC",
            "## Step6 - Within-mouse concordance",
            "## Step7 - Diagnostic Cell-style figures",
            "## Step8 - Final QC decision and next actions",
            "## Outputs and reload checks",
        ]
        for heading in required_headings:
            self.assertIn(heading, text)
        for token in [
            "PROVISIONAL_PRIMARY_FEATURES",
            "CONSERVATIVE_NO_SIGNAL_DETECTED",
            "TECHNICAL_RISK_SENSITIVITY_ONLY",
            "RAW_COMPLETE_PANEL",
            "cell_downstream_masks.tsv.gz",
            "section_downstream_decision.tsv",
            "gene_downstream_decision.tsv",
            "eos_gene_decision_summary.tsv",
            "hotspot_sensitivity_decision.tsv",
            "evidence_only_qc_release.tsv",
            "downstream_input_manifest.tsv",
            "Uncertain",
        ]:
            self.assertIn(token, text)
        self.assertNotIn("Obtain 10x poor-cycle diagnostics", text)
        self.assertNotIn("build_report.R", text)
        self.assertNotIn("report.html", text)
        self.assertNotIn("build_cell_downstream_masks(", code)
        self.assertNotRegex(code, r"(?m)^primary_include\s*<-")
        self.assertIn("does not redefine", text)
        self.assertTrue(notebook["metadata"]["kernelspec"]["name"].startswith("ir"))

    def test_region3_eos_neighbour_block_is_marked_invalid_at_point_of_use(self):
        notebook = read_notebook(NOTEBOOKS / "B1_Region3_primary_479.ipynb")
        cell_text = ["".join(cell.get("source", [])) for cell in notebook["cells"]]
        warning_indexes = [
            index
            for index, text in enumerate(cell_text)
            if "INVALID LEGACY EOS NEAREST-CELL ANALYSIS" in text
        ]
        self.assertEqual(len(warning_indexes), 1)
        warning_index = warning_indexes[0]
        warning = cell_text[warning_index]
        for token in [
            "Final_CellType_subtype_refined",
            "Final_CellType_subtype",
            "zero-distance self-matches",
            "not presently valid",
        ]:
            self.assertIn(token, warning)

        reference_index = next(
            index for index, text in enumerate(cell_text) if "reference_cells <- coords" in text
        )
        knn_index = next(
            index for index, text in enumerate(cell_text) if "nn1 <- FNN::get.knnx" in text
        )
        self.assertLess(warning_index, reference_index)
        self.assertLess(reference_index, knn_index)

    def test_region_notebook_cli_resolves_only_valid_committed_copies(self):
        renderer = REPO / "scripts" / "render_notebooks.py"
        valid = subprocess.run(
            [sys.executable, str(renderer), "--region-notebook", "Region_3"],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(valid.returncode, 0, valid.stderr)
        self.assertEqual(
            Path(valid.stdout.strip()).name,
            "01_QC_Region3.ipynb",
        )

        invalid = subprocess.run(
            [sys.executable, str(renderer), "--region-notebook", "Region_5"],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(invalid.returncode, 0)
        self.assertIn("Unknown region", invalid.stderr)

    def test_duplicate_acceptance_report_package_is_retired(self):
        report_dir = REPO / "reports" / "2026-08-15_full_hpc_qc_acceptance"
        retired_files = [
            "build_report.R",
            "report.html",
            "artifact.json",
            "candidate_cycle_affected_genes_affected_only.tsv",
        ]
        remaining = [name for name in retired_files if (report_dir / name).exists()]
        self.assertEqual(remaining, [], f"duplicate report files remain: {remaining}")


if __name__ == "__main__":
    unittest.main()
