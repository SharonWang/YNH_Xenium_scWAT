import json
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


class NotebookContracts(unittest.TestCase):
    def test_four_region_notebooks_are_parameter_locked(self):
        expected = {f"Region_{index}" for index in range(1, 5)}
        observed = set()

        for index in range(1, 5):
            path = NOTEBOOKS / f"01_section_phase0_2_QC_Region{index}.ipynb"
            self.assertTrue(path.exists(), f"missing committed notebook: {path.name}")
            notebook = read_notebook(path)
            match = re.search(
                r'^REGION_ID\s*<-\s*"(Region_[1-4])"$',
                parameter_source(notebook),
                re.MULTILINE,
            )
            self.assertIsNotNone(match, f"fixed REGION_ID missing from {path.name}")
            observed.add(match.group(1))
            self.assertEqual(notebook["metadata"]["kernelspec"]["name"], "ir")
            text = "\n".join(
                "".join(cell.get("source", [])) for cell in notebook["cells"]
            )
            for token in [
                "build_cell_downstream_masks",
                "cell_downstream_masks.tsv.gz",
                "section_downstream_decision.tsv",
                "raw objects are not modified",
            ]:
                self.assertIn(token, text)

        self.assertEqual(observed, expected)

    def test_summary_notebook_is_the_complete_reader_facing_report(self):
        notebook = read_notebook(NOTEBOOKS / "02_slide_QC_summary.ipynb")
        text = "\n".join(
            "".join(cell.get("source", [])) for cell in notebook["cells"]
        )
        required_headings = [
            "## TL;DR and QC decision",
            "## Inputs and validation",
            "## Core QC distributions",
            "## Alarm evidence and evidence-only gene tiers",
            "## Subset versus full-data burden",
            "## Spatial QC",
            "## Within-mouse concordance",
            "## Final QC decision and next actions",
        ]
        for heading in required_headings:
            self.assertIn(heading, text)
        for token in [
            "PROVISIONAL_PRIMARY_FEATURES",
            "CONSERVATIVE_NO_SIGNAL_DETECTED",
            "TECHNICAL_RISK_SENSITIVITY_ONLY",
            "RAW_COMPLETE_PANEL",
            "PENDING_DOWNSTREAM_ANALYSIS",
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
        self.assertEqual(notebook["metadata"]["kernelspec"]["name"], "ir")

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
            "01_section_phase0_2_QC_Region3.ipynb",
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
