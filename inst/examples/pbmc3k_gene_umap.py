"""Calculate group-guided gene UMAP coordinates for the PBMC3K DE result."""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import umap


def main(input_path: str, output_path: str) -> None:
    data = pd.read_csv(input_path, sep="\t")
    required = {"gene", "avg_log2FC", "p_val_adj"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing columns: {', '.join(sorted(missing))}")

    valid = (
        np.isfinite(data["avg_log2FC"])
        & np.isfinite(data["p_val_adj"])
        & (data["p_val_adj"] > 0)
    )
    data = data.loc[valid].copy()
    data["minus_log10_padj"] = -np.log10(data["p_val_adj"])
    data["group"] = "Not significant"
    data.loc[
        (data["p_val_adj"] < 0.05) & (data["avg_log2FC"] > 1), "group"
    ] = "Up"
    data.loc[
        (data["p_val_adj"] < 0.05) & (data["avg_log2FC"] < -1), "group"
    ] = "Down"

    features = np.column_stack(
        [
            data["avg_log2FC"],
            np.abs(data["avg_log2FC"]),
            data["minus_log10_padj"],
            np.sign(data["avg_log2FC"]) * data["minus_log10_padj"],
        ]
    )
    features = StandardScaler().fit_transform(features)
    group_codes = pd.Categorical(
        data["group"], categories=["Down", "Not significant", "Up"]
    ).codes

    reducer = umap.UMAP(
        n_neighbors=40,
        min_dist=0.16,
        metric="euclidean",
        target_metric="categorical",
        target_weight=0.68,
        random_state=20260825,
        n_jobs=1,
    )
    raw = reducer.fit_transform(features, y=group_codes)
    raw = StandardScaler().fit_transform(raw)

    anchors = {"Down": -2.7, "Not significant": 0.0, "Up": 2.7}
    candy_strength = 0.88
    final = np.empty_like(raw)
    for group, anchor in anchors.items():
        mask = data["group"].to_numpy() == group
        local = raw[mask] - raw[mask].mean(axis=0)
        target_x = anchor + local[:, 0] * 0.72
        y_scale = 0.72 if group == "Not significant" else 1.0
        target_y = local[:, 1] * y_scale
        final[mask, 0] = (1 - candy_strength) * raw[mask, 0] + candy_strength * target_x
        final[mask, 1] = (1 - candy_strength) * raw[mask, 1] + candy_strength * target_y

    data["UMAP_raw_1"] = raw[:, 0]
    data["UMAP_raw_2"] = raw[:, 1]
    data["UMAP_1"] = final[:, 0]
    data["UMAP_2"] = final[:, 1]
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(output_path, sep="\t", index=False)

    counts = data["group"].value_counts()
    print(f"valid_genes={len(data)}")
    for group in ["Up", "Down", "Not significant"]:
        print(f"{group.replace(' ', '_')}={int(counts.get(group, 0))}")
    print(f"output={Path(output_path).resolve()}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: pbmc3k_gene_umap.py INPUT.tsv OUTPUT.tsv")
    main(sys.argv[1], sys.argv[2])
