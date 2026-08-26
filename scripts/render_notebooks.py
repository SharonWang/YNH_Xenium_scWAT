#!/usr/bin/env python3
"""Resolve, rebuild, inject environment defaults, and validate QC notebooks."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
NOTEBOOKS = REPO / "notebooks"
REGION_NOTEBOOKS = {
    f"Region_{index}": NOTEBOOKS / f"01_QC_Region{index}.ipynb"
    for index in range(1, 5)
}


def read_notebook(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def notebook_text(notebook: dict) -> str:
    return "\n".join("".join(cell.get("source", [])) for cell in notebook["cells"])


def validate_notebook(path: Path, notebook_type: str) -> None:
    notebook = read_notebook(path)
    expected_cells = 12 if notebook_type == "region" else 10
    if len(notebook.get("cells", [])) != expected_cells:
        raise ValueError(f"{notebook_type} notebook requires {expected_cells} cells")
    expected_types = (
        ["markdown", "markdown", "code", "markdown", "code", "markdown", "code", "markdown", "code", "markdown", "code", "markdown"]
        if notebook_type == "region"
        else ["markdown", "markdown", "code", "markdown", "code", "markdown", "code", "markdown", "code", "markdown"]
    )
    if [cell.get("cell_type") for cell in notebook["cells"]] != expected_types:
        raise ValueError("notebook cell order does not match the Colon reader flow")
    ids = [cell.get("id") for cell in notebook["cells"]]
    if any(not value for value in ids) or len(ids) != len(set(ids)):
        raise ValueError("notebook cells require unique non-empty IDs")
    if notebook.get("metadata", {}).get("kernelspec", {}).get("name") != "ir":
        raise ValueError("an R kernelspec is required")
    for cell in notebook["cells"]:
        if cell.get("cell_type") == "code" and (
            cell.get("execution_count") is not None or cell.get("outputs", [])
        ):
            raise ValueError("committed QC notebooks must be output-free")
    text = notebook_text(notebook)
    if notebook_type == "region":
        required = [
            "## Goal", "## Setup", "### Key assumptions", "## Inputs and integrity",
            "## Cell QC", "## Outputs and checks", "## Next steps",
            "expected_scwat_regions()", "write_scwat_region_qc_bundle",
            "validate_scwat_region_qc_bundle", "SCWAT_QC_MODE",
        ]
        region_matches = re.findall(r'REGION_ID\s*<-\s*"(Region_[1-4])"', text)
        if len(region_matches) != 1:
            raise ValueError("region notebook requires exactly one locked REGION_ID")
    else:
        required = [
            "## Goal", "## Setup", "## Inputs and integrity", "## Cell QC",
            "## Outputs and checks", "## Next steps", "expected_scwat_regions()",
            "read_scwat_slide_qc_outputs", "summarise_scwat_mouse_sections",
            "write_scwat_slide_qc_bundle", "validate_scwat_slide_qc_bundle",
        ]
    missing = [token for token in required if token not in text]
    if missing:
        raise ValueError("missing required notebook content: " + "; ".join(missing))


def inject_environment(input_path: Path, output_path: Path, assignments: list[str]) -> None:
    """Prepend `Sys.setenv()` overrides without editing the locked notebook body."""
    notebook = read_notebook(input_path)
    pairs = []
    for assignment in assignments:
        if "=" not in assignment:
            raise ValueError(f"invalid --set value: {assignment}")
        key, value = assignment.split("=", 1)
        if not re.fullmatch(r"SCWAT_[A-Z0-9_]+", key):
            raise ValueError(f"only SCWAT_* environment variables may be injected: {key}")
        escaped = value.replace("\\", "\\\\").replace('"', '\\"')
        pairs.append(f'{key} = "{escaped}"')
    source = "Sys.setenv(\n  " + ",\n  ".join(pairs) + "\n)\n"
    notebook["cells"].insert(0, {
        "cell_type": "code", "execution_count": None, "id": "runtime-environment",
        "metadata": {"tags": ["injected-runtime"]}, "outputs": [], "source": source.splitlines(keepends=True),
    })
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(notebook, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--region-notebook", choices=sorted(REGION_NOTEBOOKS))
    parser.add_argument("--validate", type=Path)
    parser.add_argument("--type", choices=["region", "summary"], default="region")
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--inject", nargs=2, metavar=("INPUT", "OUTPUT"))
    parser.add_argument("--set", action="append", default=[])
    args = parser.parse_args()

    if args.build:
        subprocess.run([sys.executable, str(REPO / "scripts" / "build_scwat_initial_qc_notebooks.py")], check=True)
        return
    if args.region_notebook:
        path = REGION_NOTEBOOKS[args.region_notebook]
        validate_notebook(path, "region")
        print(path.resolve())
        return
    if args.validate:
        validate_notebook(args.validate, args.type)
        print(f"Notebook validation passed: {args.validate}")
        return
    if args.inject:
        inject_environment(Path(args.inject[0]), Path(args.inject[1]), args.set)
        return
    parser.error("select --build, --region-notebook, --validate, or --inject")


if __name__ == "__main__":
    main()
