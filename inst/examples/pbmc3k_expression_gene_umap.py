"""Unsupervised gene UMAP from a gene-by-cell expression matrix."""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.io import mmread
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap


def main(
    matrix_path: str,
    gene_path: str,
    de_path: str,
    output_path: str,
    n_components: int = 2,
) -> None:
    matrix = mmread(matrix_path).tocsr()
    genes = pd.read_csv(gene_path, sep="\t")
    de = pd.read_csv(de_path, sep="\t")
    if matrix.shape[0] != len(genes):
        raise ValueError("Gene-name count does not match expression-matrix rows.")

    # A gene is one observation; its expression across 800 cells is its feature
    # vector. Row-wise z-scoring emphasizes expression pattern over abundance.
    expression = matrix.toarray().astype(np.float32, copy=False)
    row_mean = expression.mean(axis=1, keepdims=True)
    row_sd = expression.std(axis=1, keepdims=True)
    informative = row_sd[:, 0] > 0
    expression = (expression[informative] - row_mean[informative]) / row_sd[informative]
    retained_genes = genes.loc[informative, "gene"].astype(str).reset_index(drop=True)

    pca = PCA(n_components=30, svd_solver="randomized", random_state=20260825)
    principal_components = pca.fit_transform(expression)
    principal_components = StandardScaler().fit_transform(principal_components)

    reducer = umap.UMAP(
        n_neighbors=35,
        min_dist=0.22,
        metric="cosine",
        n_components=n_components,
        random_state=20260825,
        n_jobs=1,
    )
    embedding = reducer.fit_transform(principal_components)

    coordinates = pd.DataFrame({"gene": retained_genes})
    for component in range(n_components):
        coordinates[f"UMAP_{component + 1}"] = embedding[:, component]
    de_columns = de[["gene", "avg_log2FC", "p_val_adj"]].drop_duplicates("gene")
    result = coordinates.merge(de_columns, on="gene", how="inner", validate="one_to_one")
    valid = (
        np.isfinite(result["avg_log2FC"])
        & np.isfinite(result["p_val_adj"])
        & (result["p_val_adj"] > 0)
    )
    result = result.loc[valid].copy()
    result["group"] = "Not significant"
    result.loc[
        (result["p_val_adj"] < 0.05) & (result["avg_log2FC"] > 1), "group"
    ] = "Up"
    result.loc[
        (result["p_val_adj"] < 0.05) & (result["avg_log2FC"] < -1), "group"
    ] = "Down"
    result.to_csv(output_path, sep="\t", index=False)

    counts = result["group"].value_counts()
    print(f"matrix_genes={matrix.shape[0]}")
    print(f"matrix_cells={matrix.shape[1]}")
    print(f"genes_with_variable_expression={informative.sum()}")
    print(f"umap_dimensions={n_components}")
    print(f"pca_variance_30={pca.explained_variance_ratio_.sum():.6f}")
    for group in ["Up", "Down", "Not significant"]:
        print(f"{group.replace(' ', '_')}={int(counts.get(group, 0))}")
    print(f"output={Path(output_path).resolve()}")


if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        raise SystemExit(
            "Usage: pbmc3k_expression_gene_umap.py MATRIX.mtx GENES.tsv "
            "DE.tsv OUTPUT.tsv [N_COMPONENTS]"
        )
    dimensions = int(sys.argv[5]) if len(sys.argv) == 6 else 2
    if dimensions not in (2, 3):
        raise SystemExit("N_COMPONENTS must be 2 or 3.")
    main(*sys.argv[1:5], n_components=dimensions)
