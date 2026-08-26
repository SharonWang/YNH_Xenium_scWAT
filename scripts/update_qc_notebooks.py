"""Apply the approved fixed-QC and Eos-warning edits to committed notebooks."""

from __future__ import annotations

import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
NOTEBOOKS = REPO / "notebooks"


def read_notebook(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_notebook(path: Path, notebook: dict) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(notebook, handle, indent=1, ensure_ascii=False)
        handle.write("\n")


def source_text(cell: dict) -> str:
    source = cell.get("source", [])
    return source if isinstance(source, str) else "".join(source)


def set_source(cell: dict, text: str) -> None:
    cell["source"] = text.splitlines(keepends=True)


def replace_once(text: str, old: str, new: str, context: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"Expected one {context} match, found {count}")
    return text.replace(old, new, 1)


def edit_region_notebook(path: Path) -> dict:
    notebook = read_notebook(path)
    parameter_cells = [
        cell
        for cell in notebook["cells"]
        if "parameters" in cell.get("metadata", {}).get("tags", [])
    ]
    if len(parameter_cells) != 1:
        raise RuntimeError(f"{path.name}: expected one parameter cell")

    parameter_cell = parameter_cells[0]
    parameters = source_text(parameter_cell)
    if "FIXED_CELL_QC_THRESHOLDS_PATH" not in parameters:
        parameters = (
            parameters.rstrip()
            + '\nFIXED_CELL_QC_THRESHOLDS_PATH <- file.path(PIPELINE_REPO, "config", "fixed_cell_qc_thresholds.tsv")\n'
        )
    set_source(parameter_cell, parameters)

    found_source_call = False
    found_qc_call = False
    found_qc_explanation = False
    found_mask_explanation = False
    found_transcript_initialization = False
    for cell in notebook["cells"]:
        text = source_text(cell)
        if cell.get("cell_type") == "code" and 'source(file.path(PIPELINE_REPO, "R", "source.R"))' in text:
            if "read_fixed_cell_qc_thresholds" not in text:
                text = text.rstrip() + "\nfixed_thresholds <- read_fixed_cell_qc_thresholds(FIXED_CELL_QC_THRESHOLDS_PATH)\n"
            found_source_call = True
        if (
            cell.get("cell_type") == "code"
            and 'if (extended_mode == "FULL_HPC") {' in text
            and "summarise_transcript_quality_arrow" in text
        ):
            if "transcript_gene_quality <- data.frame()" not in text:
                text = "transcript_gene_quality <- data.frame()\n" + text
            found_transcript_initialization = True
        if "qc <- calculate_xenium_cell_qc(xenium$counts, xenium$cells, REGION_ID)" in text:
            text = text.replace(
                "qc <- calculate_xenium_cell_qc(xenium$counts, xenium$cells, REGION_ID)",
                "qc <- calculate_xenium_cell_qc(\n"
                "  xenium$counts, xenium$cells, REGION_ID,\n"
                "  fixed_thresholds = fixed_thresholds\n"
                ")",
            )
            found_qc_call = True
        elif "fixed_thresholds = fixed_thresholds" in text and "calculate_xenium_cell_qc" in text:
            found_qc_call = True
        if cell.get("cell_type") == "markdown" and "qc_core_pass = TRUE" in text:
            text = (
                "### Prespecified core cell-QC rule\n\n"
                "A cell passes `qc_core_pass` only when both strict open intervals hold:\n\n"
                "- `5 < nFeature_Xenium < 200`\n"
                "- `10 < nCount_Xenium < 1000`\n\n"
                "Values exactly equal to 5, 200, 10, or 1000 fail. Robust distribution limits remain diagnostic only and do not redefine this cohort.\n"
            )
            found_qc_explanation = True
        elif cell.get("cell_type") == "markdown" and "Prespecified core cell-QC rule" in text:
            found_qc_explanation = True
        if cell.get("cell_type") == "markdown" and "primary_include: passes the core QC" in text:
            text = (
                "### Downstream masks (raw objects remain unchanged)\n\n"
                "- `primary_include`: exactly `qc_core_pass`, using `5 < nFeature_Xenium < 200` and `10 < nCount_Xenium < 1000`.\n"
                "- `strict_include`: `primary_include & !qc_review_flag`; this additionally removes cells carrying any nucleus, segmentation, area, control, or core-QC review flag.\n"
                "- `hotspot_sensitivity_include`: `primary_include` with Region 3 morphology-review hotspot cells excluded.\n\n"
                "The primary mask is the main downstream cohort. All flags and all raw cells remain stored for audit and sensitivity analyses.\n"
            )
            found_mask_explanation = True
        elif cell.get("cell_type") == "markdown" and "Downstream masks (raw objects remain unchanged)" in text:
            found_mask_explanation = True
        set_source(cell, text)

    if not all((found_source_call, found_qc_call, found_qc_explanation, found_mask_explanation, found_transcript_initialization)):
        raise RuntimeError(
            f"{path.name}: missing target(s): source_call={found_source_call}, qc_call={found_qc_call}, "
            f"qc_explanation={found_qc_explanation}, mask_explanation={found_mask_explanation}, "
            f"transcript_initialization={found_transcript_initialization}"
        )
    return notebook


def edit_summary_notebook(path: Path) -> dict:
    notebook = read_notebook(path)
    setup_updated = False
    decision_updated = False
    validation_updated = False
    for cell in notebook["cells"]:
        text = source_text(cell)
        if cell.get("cell_type") == "code" and "SUBSET_REFERENCE_PATH <-" in text:
            if "FIXED_CELL_QC_THRESHOLDS_PATH" not in text:
                text = (
                    text.rstrip()
                    + '\nFIXED_CELL_QC_THRESHOLDS_PATH <- file.path(PIPELINE_REPO, "config", "fixed_cell_qc_thresholds.tsv")\n'
                )
            setup_updated = True
        if cell.get("cell_type") == "code" and 'source(file.path(PIPELINE_REPO, "R", "source.R"))' in text:
            if "read_fixed_cell_qc_thresholds" not in text:
                text = text.rstrip() + "\nfixed_thresholds <- read_fixed_cell_qc_thresholds(FIXED_CELL_QC_THRESHOLDS_PATH)\n"
        if cell.get("cell_type") == "markdown" and text.strip() == "## Step2 - QC decision":
            text = (
                "## Step2 - QC decision\n\n"
                "`primary_include` is imported from each region artifact and does not redefine the mask at slide level. "
                "The notebook verifies that every imported value equals the prespecified strict rule "
                "`5 < nFeature_Xenium < 200` and `10 < nCount_Xenium < 1000`.\n"
            )
            decision_updated = True
        elif cell.get("cell_type") == "markdown" and "does not redefine the mask at slide level" in text:
            decision_updated = True
        if cell.get("cell_type") == "markdown" and text.strip() == "## Step8 - Final QC decision and next actions":
            text = (
                "## Step8 - Final QC decision and next actions\n\n"
                "Region 4 remains mapping-only and cannot influence PCA, integration, clustering, marker selection, or reference labels. "
                "Its later mapped cells must retain `Uncertain` whenever assignment confidence is insufficient.\n"
            )
        validation_block = (
            "expected_primary_include <- with(\n"
            "  evidence_summary$cell_masks,\n"
            "  apply_fixed_primary_bounds(nFeature_Xenium, nCount_Xenium, fixed_thresholds)\n"
            ")\n"
            "stopifnot(identical(\n"
            "  as.logical(evidence_summary$cell_masks$primary_include),\n"
            "  as.logical(expected_primary_include)\n"
            "))"
        )
        if cell.get("cell_type") == "code" and "dim(evidence_summary$cell_masks)" in text:
            while validation_block in text:
                text = text.replace("\n" + validation_block, "", 1)
            text = text.rstrip() + "\n" + validation_block + "\n"
            validation_updated = True
        set_source(cell, text)
    if not all((setup_updated, decision_updated, validation_updated)):
        raise RuntimeError(
            f"{path.name}: incomplete edit: setup={setup_updated}, "
            f"decision={decision_updated}, validation={validation_updated}"
        )
    return notebook


def edit_eos_notebook(path: Path) -> dict:
    notebook = read_notebook(path)
    if any("INVALID LEGACY EOS NEAREST-CELL ANALYSIS" in source_text(cell) for cell in notebook["cells"]):
        return notebook
    reference_index = next(
        index
        for index, cell in enumerate(notebook["cells"])
        if "reference_cells <- coords" in source_text(cell)
    )
    warning = {
        "cell_type": "markdown",
        "id": "invalid-eos-nearest-warning",
        "metadata": {},
        "source": [],
    }
    set_source(
        warning,
        "### ⚠ INVALID LEGACY EOS NEAREST-CELL ANALYSIS\n\n"
        "The following spatial-neighbour block is retained only to document a known error and its existing outputs are **not presently valid**. "
        "Rescued Eos cells were assigned with `Final_CellType_subtype_refined`, but `celltype_col` and the non-Eos `reference_cells` pool use the older `Final_CellType_subtype`. "
        "Some rescued Eos therefore occur in both `eos_all` and `reference_cells`, allowing `FNN::get.knnx()` to return zero-distance self-matches. "
        "Do not interpret `closest_df`, `nearest_by_type`, or neighbourhood summaries until query and reference cell IDs are explicitly disjoint and both are defined from the same frozen refined label.\n",
    )
    notebook["cells"].insert(reference_index, warning)

    reference_cell = notebook["cells"][reference_index + 1]
    set_source(
        reference_cell,
        "# KNOWN INVALID LEGACY LOGIC: celltype_col points to the unrefined label.\n"
        "# Rescued Eos can therefore remain in this apparent non-Eos pool.\n"
        + source_text(reference_cell),
    )
    knn_cell = next(
        cell for cell in notebook["cells"] if "nn1 <- FNN::get.knnx" in source_text(cell)
    )
    set_source(
        knn_cell,
        "# NOT VALID FOR INTERPRETATION: eos_all and reference_cells were not\n"
        "# checked for disjoint cell IDs, so zero-distance self-matches occur.\n"
        + source_text(knn_cell),
    )
    return notebook


def main() -> None:
    edited = {}
    for region in range(1, 5):
        path = NOTEBOOKS / f"01_QC_Region{region}.ipynb"
        edited[path] = edit_region_notebook(path)
    summary_path = NOTEBOOKS / "02_slide_QC_summary.ipynb"
    eos_path = NOTEBOOKS / "B1_Region3_primary_479.ipynb"
    edited[summary_path] = edit_summary_notebook(summary_path)
    edited[eos_path] = edit_eos_notebook(eos_path)
    for path, notebook in edited.items():
        write_notebook(path, notebook)
    print("Updated four region notebooks, slide summary, and Region 3 Eos warning.")


if __name__ == "__main__":
    main()
