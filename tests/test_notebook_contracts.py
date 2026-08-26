import json
import subprocess
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
NOTEBOOKS = REPO / "notebooks"

def read_notebook(path):
    return json.loads(path.read_text(encoding="utf-8"))

def notebook_text(notebook):
    return "\n".join("".join(cell.get("source", [])) for cell in notebook["cells"])

class NotebookContracts(unittest.TestCase):
    def test_initial_qc_notebooks_match_the_colon_reader_flow(self):
        expected_types = ["markdown", "markdown", "code", "markdown", "code", "markdown", "code", "markdown", "code", "markdown", "code", "markdown"]
        for index in range(1, 5):
            notebook = read_notebook(NOTEBOOKS / f"01_QC_Region{index}.ipynb")
            self.assertEqual(len(notebook["cells"]), 12)
            self.assertEqual([cell["cell_type"] for cell in notebook["cells"]], expected_types)
            text = notebook_text(notebook)
            self.assertIn(f'REGION_ID <- "Region_{index}"', text)
            for token in ["## Goal", "## Setup", "### Key assumptions", "## Inputs and integrity", "## Cell QC", "## Outputs and checks", "## Next steps", "SCWAT_QC_MODE", "scwat_sample_manifest.tsv", "expected_scwat_regions()", "write_scwat_region_qc_bundle", "validate_scwat_region_qc_bundle", "adipose_data", "adipose_analysis"]:
                self.assertIn(token, text)
            self.assertNotIn("COLON_", text)
            self.assertNotIn("colon_analysis", text)
            self.assertEqual(notebook["metadata"]["kernelspec"]["name"], "ir")
            self.assertTrue(all(cell.get("execution_count") is None and cell.get("outputs", []) == [] for cell in notebook["cells"] if cell["cell_type"] == "code"))

    def test_slide_summary_matches_the_colon_reader_flow_for_four_scwat_regions(self):
        notebook = read_notebook(NOTEBOOKS / "02_slide_QC_summary.ipynb")
        self.assertEqual(len(notebook["cells"]), 10)
        self.assertEqual([cell["cell_type"] for cell in notebook["cells"]], ["markdown", "markdown", "code", "markdown", "code", "markdown", "code", "markdown", "code", "markdown"])
        text = notebook_text(notebook)
        for token in ["expected_scwat_regions()", "read_scwat_slide_qc_outputs", "summarise_scwat_mouse_sections", "write_scwat_slide_qc_bundle", "validate_scwat_slide_qc_bundle", "Mouse_1", "Mouse_2"]:
            self.assertIn(token, text)
        self.assertNotIn("position", text.lower())
        self.assertNotIn("manifest$side", text)
        self.assertIn("without inventing a left/right factor", text)
        self.assertTrue(all(cell.get("execution_count") is None and cell.get("outputs", []) == [] for cell in notebook["cells"] if cell["cell_type"] == "code"))

    def test_region3_eos_neighbour_block_remains_marked_invalid(self):
        notebook = read_notebook(NOTEBOOKS / "B1_Region3_primary_479.ipynb")
        cells = ["".join(cell.get("source", [])) for cell in notebook["cells"]]
        warning_index = next(i for i, text in enumerate(cells) if "INVALID LEGACY EOS NEAREST-CELL ANALYSIS" in text)
        for token in ["Final_CellType_subtype_refined", "Final_CellType_subtype", "zero-distance self-matches", "not presently valid"]:
            self.assertIn(token, cells[warning_index])
        reference_index = next(i for i, text in enumerate(cells) if "reference_cells <- coords" in text)
        knn_index = next(i for i, text in enumerate(cells) if "nn1 <- FNN::get.knnx" in text)
        self.assertLess(warning_index, reference_index)
        self.assertLess(reference_index, knn_index)

    def test_renderer_resolves_and_validates_only_four_regions(self):
        renderer = REPO / "scripts" / "render_notebooks.py"
        for index in range(1, 5):
            result = subprocess.run([sys.executable, str(renderer), "--region-notebook", f"Region_{index}"], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
        invalid = subprocess.run([sys.executable, str(renderer), "--region-notebook", "Region_5"], capture_output=True, text=True)
        self.assertNotEqual(invalid.returncode, 0)

    def test_duplicate_acceptance_report_package_is_retired(self):
        report_dir = REPO / "reports" / "2026-08-15_full_hpc_qc_acceptance"
        retired = ["build_report.R", "report.html", "artifact.json", "candidate_cycle_affected_genes_affected_only.tsv"]
        self.assertEqual([name for name in retired if (report_dir / name).exists()], [])

if __name__ == "__main__":
    unittest.main()
