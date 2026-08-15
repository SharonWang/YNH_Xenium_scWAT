import json
import re
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
            "## Alarm evidence and candidate genes",
            "## Subset versus full-data burden",
            "## Spatial QC",
            "## Within-mouse concordance",
            "## Final QC decision and next actions",
        ]
        for heading in required_headings:
            self.assertIn(heading, text)
        self.assertIn("CANDIDATE_NOT_CONFIRMED", text)
        self.assertIn("REQUIRES_10X_DIAGNOSTICS", text)
        self.assertNotIn("build_report.R", text)
        self.assertNotIn("report.html", text)
        self.assertEqual(notebook["metadata"]["kernelspec"]["name"], "ir")


if __name__ == "__main__":
    unittest.main()
