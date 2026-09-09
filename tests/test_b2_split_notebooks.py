import json
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
NOTEBOOKS = REPO / "notebooks"
BRANCHES = ("all_QCpass", "adipose_only", "lymph_node_only")


def read_notebook(path):
    return json.loads(path.read_text(encoding="utf-8"))


def notebook_text(notebook):
    return "\n".join("".join(cell.get("source", [])) for cell in notebook["cells"])


class SplitB2NotebookContracts(unittest.TestCase):
    def test_twelve_region_branch_notebooks_are_clean_and_region_locked(self):
        for region in range(1, 5):
            for branch in BRANCHES:
                path = NOTEBOOKS / f"B2_Region{region}_{branch}_479.ipynb"
                self.assertTrue(path.exists(), path)
                notebook = read_notebook(path)
                text = notebook_text(notebook)
                self.assertIn(f'REGION_ID <- "Region_{region}"', text)
                self.assertIn(f'ANALYSIS_BRANCH <- "{branch}"', text)
                self.assertIn("primary_include_revised", text)
                self.assertIn("FIX_GENESET <- raw_panel_genes", text)
                self.assertIn("stopifnot(length(FIX_GENESET) == 479L)", text)
                self.assertIn("dims_use <- 1:30", text)
                self.assertIn("LogNormalize", text)
                self.assertIn("summarise_pca_qc_correlations", text)
                self.assertIn("summarise_cluster_stability", text)
                self.assertIn("CLUSTER_ALGORITHM", text)
                self.assertNotIn("algorithm = 4,", text)
                self.assertIn("Uncertain", text)
                self.assertIn("build_eos_spatial_pools", text)
                self.assertIn("calculate_eos_knn_edges", text)
                self.assertIn("optional_packages", text)
                self.assertIn("RUN_DIPTEST", text)
                self.assertIn("RUN_MCLUST", text)
                self.assertIn("RUN_CELLCHAT", text)
                self.assertIn("RUN_WANG_INTEGRATION", text)
                self.assertIn("run_mclust_diagnostic", text)
                self.assertNotIn("mclust::Mclust(", text)
                self.assertIn("apply_scwat_cell_type_order", text)
                self.assertIn("order_marker_features", text)
                self.assertIn('group.by = "Final_CellType_subtype"', text)
                self.assertIn("## 10.1", text)
                self.assertIn("## 10.2", text)
                self.assertIn("## 10.3", text)
                self.assertIn("## 10.4", text)
                self.assertIn("build_eos_spatial_pools", text)
                self.assertIn("calculate_eos_knn_edges", text)
                self.assertIn("calculate_eos_distance_by_cell_type", text)
                self.assertIn("rank_eos_state_knn_associations", text)
                self.assertIn("run_eos_spatial_cellchat", text)
                self.assertIn("filter_eos_cellchat_interactions", text)
                self.assertIn("run_wang_xenium_integration", text)
                self.assertIn("summarise_cross_dataset_eos_neighbours", text)
                self.assertIn("CELLCHAT_STATE_PROPORTION <- 0.30", text)
                self.assertIn("cellchat_significant_interactions", text)
                self.assertIn("wang_xenium_integration_status", text)
                self.assertIn("wang_xenium_shared_gene_manifest.tsv", text)
                self.assertIn("wang_xenium_cross_dataset_neighbours.tsv", text)
                self.assertIn("eos_knn_k1_edges.tsv", text)
                self.assertIn("eos_knn_k15_edges.tsv.gz", text)
                self.assertIn("eos_knn_state_association.tsv", text)
                self.assertIn("eos_spatial_cellchat_status.tsv", text)
                self.assertIn("eos_spatial_cellchat_significant.tsv", text)
                self.assertNotIn("EXPLORATORY_SPATIAL_COEXPRESSION_NOT_CELLCHAT_INFERENCE", text)
                self.assertIn("write_validated_seurat_checkpoint", text)
                self.assertLess(
                    text.index("validate_runtime_paths"),
                    text.index("dir.create(REGION_SHARED_ROOT")
                )
                self.assertNotIn("install.packages", text)
                self.assertNotIn("\ufffd", text)
                self.assertNotIn("C:\\\\", text)
                self.assertEqual(notebook["metadata"]["kernelspec"]["name"], "ir")
                self.assertTrue(all(
                    cell.get("execution_count") is None and cell.get("outputs", []) == []
                    for cell in notebook["cells"] if cell["cell_type"] == "code"
                ))

    def test_all_qcpass_notebooks_freeze_the_partition_for_child_branches(self):
        for region in range(1, 5):
            notebook = read_notebook(NOTEBOOKS / f"B2_Region{region}_all_QCpass_479.ipynb")
            text = notebook_text(notebook)
            self.assertIn("build_tissue_branch_manifest", text)
            self.assertIn("lymph_node_domain_manifest.tsv.gz", text)
            self.assertIn("lymph_node_boundary_sensitivity.tsv", text)

    def test_child_notebooks_consume_but_do_not_rederive_the_frozen_partition(self):
        for region in range(1, 5):
            for branch in ("adipose_only", "lymph_node_only"):
                notebook = read_notebook(NOTEBOOKS / f"B2_Region{region}_{branch}_479.ipynb")
                text = notebook_text(notebook)
                self.assertIn("lymph_node_domain_manifest.tsv.gz", text)
                self.assertIn("select_tissue_branch_ids", text)
                self.assertNotIn("derive_lymph_node_domain(", text)

    def test_region4_is_visibly_sensitivity_only(self):
        for branch in BRANCHES:
            notebook = read_notebook(NOTEBOOKS / f"B2_Region4_{branch}_479.ipynb")
            text = notebook_text(notebook)
            self.assertIn('SECTION_ROLE <- "SENSITIVITY_ONLY"', text)
            self.assertIn("must not define the final reference", text)


if __name__ == "__main__":
    unittest.main()
