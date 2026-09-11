#!/usr/bin/env python3
"""
ATLAS PBMC RNA-ATAC preprocessing pipeline.

Converted from the original PBMC data.ipynb while retaining the complete
raw-data preprocessing and ratio-split workflow:

1. read the 10x Multiome filtered feature-barcode H5;
2. split Gene Expression and Peaks;
3. preprocess RNA and write RNA_counts_qc.h5ad;
4. preprocess ATAC and write ATAC_counts_qc.h5ad;
5. construct ATAC gene activity from the GTF and write ATAC_gas.h5ad;
6. construct feature_aligned.h5ad;
7. generate train/validation split files for multiple single-modality ratios.

The three core intermediate files are regenerated from the raw 10x H5 on every run.
The fragments TSV.GZ path is retained in the config for provenance / future use,
but the original notebook preprocessing code does not consume it.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.model_selection import train_test_split

from common import (
    ensure_dir,
    get_common_names,
    get_project_root,
    load_config,
    require_file,
    require_obs_column,
    resolve_project_path,
    safe_write_h5ad,
    subset_and_copy,
)


def set_seed(seed: int = 1234) -> None:
    random.seed(seed)
    np.random.seed(seed)


def tfidf_transform_inplace(adata: ad.AnnData, scale_factor: float = 1e4) -> None:
    """Retain the custom TF-IDF transform used in the original PBMC notebook."""
    if sparse.issparse(adata.X):
        X = adata.X.tocsr().astype(np.float32)
    else:
        X = sparse.csr_matrix(np.asarray(adata.X, dtype=np.float32))

    cell_sums = np.asarray(X.sum(axis=1)).ravel()
    cell_sums[cell_sums == 0] = 1.0
    tf = X.multiply(1.0 / cell_sums[:, None])

    n_cells = X.shape[0]
    feature_counts = np.asarray((X > 0).sum(axis=0)).ravel()
    idf = np.log1p(n_cells / (1.0 + feature_counts)).astype(np.float32)

    X_tfidf = tf.multiply(idf)
    X_tfidf = X_tfidf * scale_factor
    X_tfidf.data = np.log1p(X_tfidf.data)
    adata.X = X_tfidf.tocsr()


def normalize_log1p_inplace(adata: ad.AnnData, target_sum: float = 1e4) -> None:
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)


def add_counts_layer(adata: ad.AnnData) -> None:
    adata.layers["counts"] = adata.X.copy()


def choose_hvg_inplace(
    adata: ad.AnnData,
    batch_key: str | None = None,
    min_mean: float = 0.0125,
    max_mean: float = 3.0,
    min_disp: float = 0.5,
) -> None:
    kwargs = {
        "min_mean": min_mean,
        "max_mean": max_mean,
        "min_disp": min_disp,
    }
    if batch_key is not None and batch_key in adata.obs.columns:
        kwargs["batch_key"] = batch_key
    sc.pp.highly_variable_genes(adata, **kwargs)


@dataclass
class PartialModalityAssignment:
    modality: Dict[str, str]
    paired_cells: List[str]
    rna_only_cells: List[str]
    atac_only_cells: List[str]


def assign_partial_modality(
    cells: Sequence[str],
    single_frac: float = 0.2,
    rna_keep_prob: float = 0.5,
    seed: int = 1234,
) -> PartialModalityAssignment:
    """Retain the partial-pairing assignment logic from the notebook."""
    rng = np.random.default_rng(seed)
    cells = np.asarray(list(cells), dtype=object)
    n_single = int(round(len(cells) * single_frac))

    if n_single > 0:
        single_cells = rng.choice(cells, size=n_single, replace=False)
    else:
        single_cells = np.array([], dtype=object)

    paired_cells = sorted(set(cells.tolist()) - set(single_cells.tolist()))
    modality = {c: "paired" for c in cells.tolist()}

    if n_single > 0:
        single_is_rna = rng.random(n_single) < rna_keep_prob
        for cell, keep_rna in zip(single_cells.tolist(), single_is_rna.tolist()):
            modality[cell] = "RNA" if keep_rna else "ATAC"

    rna_only = sorted([c for c, m in modality.items() if m == "RNA"])
    atac_only = sorted([c for c, m in modality.items() if m == "ATAC"])

    return PartialModalityAssignment(
        modality=modality,
        paired_cells=paired_cells,
        rna_only_cells=rna_only,
        atac_only_cells=atac_only,
    )


def read_10x_multiome_h5(input_h5: Path) -> Tuple[ad.AnnData, ad.AnnData]:
    """
    Read the 10x Multiome filtered feature-barcode H5 and split it into
    RNA/Gene Expression and ATAC/Peaks, matching the original notebook.
    """
    adata_all = sc.read_10x_h5(str(input_h5), gex_only=False)

    if "feature_types" not in adata_all.var.columns:
        raise ValueError(
            "The 10x H5 does not contain var['feature_types']; "
            "cannot split Gene Expression and Peaks."
        )

    feature_types = adata_all.var["feature_types"].astype(str)
    rna_mask = feature_types.eq("Gene Expression")
    atac_mask = feature_types.eq("Peaks")

    if int(rna_mask.sum()) == 0:
        raise ValueError("No 'Gene Expression' features found in the 10x H5.")
    if int(atac_mask.sum()) == 0:
        raise ValueError("No 'Peaks' features found in the 10x H5.")

    rna = adata_all[:, rna_mask.values].copy()
    atac = adata_all[:, atac_mask.values].copy()

    rna.var_names_make_unique()
    atac.var_names_make_unique()

    common_cells = np.intersect1d(rna.obs_names, atac.obs_names)
    rna = rna[common_cells].copy()
    atac = atac[common_cells].copy()

    return rna, atac


def preprocess_rna(
    rna: ad.AnnData,
    cfg: Dict,
) -> ad.AnnData:
    """
    PBMC RNA preprocessing exactly follows the original notebook:
    QC metrics -> raw counts layer -> normalize_total -> log1p -> HVG selection.

    Note: the original PBMC notebook calculates QC metrics but does NOT
    filter cells/genes using QC thresholds here. This behavior is preserved.
    """
    p = cfg["rna"]
    rna = rna.copy()

    batch_key = str(p["batch_key"])
    require_obs_column(rna, batch_key, "batch0")

    if "mt" not in rna.var.columns:
        rna.var["mt"] = rna.var_names.astype(str).str.startswith(str(p["mt_prefix"]))

    sc.pp.calculate_qc_metrics(
        rna,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    add_counts_layer(rna)
    normalize_log1p_inplace(rna, target_sum=float(p["target_sum"]))
    choose_hvg_inplace(
        rna,
        batch_key=batch_key,
        min_mean=float(p["hvg_min_mean"]),
        max_mean=float(p["hvg_max_mean"]),
        min_disp=float(p["hvg_min_disp"]),
    )
    return rna


def preprocess_atac(
    atac: ad.AnnData,
    cfg: Dict,
) -> ad.AnnData:
    """
    PBMC ATAC preprocessing exactly follows the original notebook:
    QC metrics -> raw counts layer -> custom TF-IDF -> HVG selection.

    Note: the original PBMC notebook does NOT perform explicit peak/cell
    filtering in this step. This behavior is preserved.
    """
    p = cfg["atac"]
    batch_key = str(cfg["rna"]["batch_key"])
    atac = atac.copy()

    require_obs_column(atac, batch_key, "batch0")
    sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)

    add_counts_layer(atac)
    tfidf_transform_inplace(atac, scale_factor=float(p["tfidf_scale_factor"]))

    choose_hvg_inplace(
        atac,
        batch_key=batch_key,
        min_mean=float(p["hvg_min_mean"]),
        max_mean=float(p["hvg_max_mean"]),
        min_disp=float(p["hvg_min_disp"]),
    )
    return atac


def build_gene_activity(
    atac_qc_path: Path,
    gtf_path: Path,
    cfg: Dict,
) -> ad.AnnData:
    """
    Construct gene-level activity from peak-level ATAC using EpiScanpy.

    The sequence is intentionally retained from the notebook:
      - save TF-IDF-transformed X into layers['normalized'];
      - reset X to raw peak counts;
      - run epi.tl.geneactivity(..., annotation='HAVANA');
      - remove duplicated gene names;
      - preserve gene-activity counts;
      - normalize_total + log1p;
      - select HVGs.
    """
    import episcanpy as epi

    p = cfg["gene_activity"]
    atac = sc.read_h5ad(str(atac_qc_path))

    if "counts" not in atac.layers:
        raise ValueError(
            "ATAC_counts_qc.h5ad does not contain layers['counts'] "
            "needed for gene activity."
        )

    atac.layers["normalized"] = atac.X.copy()
    atac.X = atac.layers["counts"].copy()

    atac_gas = epi.tl.geneactivity(
        atac,
        str(gtf_path),
        annotation=str(p["annotation"]),
    )
    if atac_gas is None:
        # Some EpiScanpy versions may mutate and return None.
        atac_gas = atac

    atac_gas = atac_gas[:, ~atac_gas.var_names.duplicated()].copy()
    require_obs_column(atac_gas, str(cfg["rna"]["batch_key"]), "batch0")

    # Preserve the unnormalized gene-activity matrix before normalization.
    add_counts_layer(atac_gas)
    normalize_log1p_inplace(
        atac_gas,
        target_sum=float(p["target_sum"]),
    )
    choose_hvg_inplace(
        atac_gas,
        batch_key=str(cfg["rna"]["batch_key"]),
        min_mean=float(p["hvg_min_mean"]),
        max_mean=float(p["hvg_max_mean"]),
        min_disp=float(p["hvg_min_disp"]),
    )

    return atac_gas


def build_feature_aligned(
    rna_qc_path: Path,
    atac_gas_path: Path,
) -> ad.AnnData:
    """
    Build the same feature-aligned object used in the original PBMC notebook.

    aligned_genes =
        union(RNA_HVG, ATAC_activity_HVG)
        intersect RNA genes
        intersect ATAC gene-activity genes

    Therefore every retained feature is present in both modalities, and is an
    HVG in at least one of the two modalities.
    """
    rna = sc.read_h5ad(str(rna_qc_path))
    atac_gas = sc.read_h5ad(str(atac_gas_path))

    rna_hvg = rna.var_names[rna.var["highly_variable"].fillna(False)].astype(str).tolist()
    atac_hvg = (
        atac_gas.var_names[atac_gas.var["highly_variable"].fillna(False)]
        .astype(str)
        .tolist()
    )

    union_hvg = set(rna_hvg) | set(atac_hvg)
    aligned_genes = sorted(
        union_hvg
        & set(rna.var_names.astype(str))
        & set(atac_gas.var_names.astype(str))
    )

    if len(aligned_genes) == 0:
        raise ValueError(
            "No common aligned genes between RNA HVGs and "
            "ATAC gene-activity HVGs."
        )

    rna_sub = rna[:, aligned_genes].copy()
    atac_sub = atac_gas[:, aligned_genes].copy()

    rna_sub.obs["modality"] = "RNA"
    atac_sub.obs["modality"] = "ATAC_gene_activity"

    aligned = ad.concat(
        [rna_sub, atac_sub],
        join="inner",
        label="modality_concat",
        keys=["RNA", "ATAC_gene_activity"],
        index_unique="__",
    )
    aligned.uns["rna_hvg"] = rna_hvg
    aligned.uns["atac_hvg"] = atac_hvg
    aligned.uns["aligned_genes"] = aligned_genes
    return aligned


def save_split_h5ads(
    outdir: Path,
    train_rna_ref: ad.AnnData,
    train_atac_full: ad.AnnData,
    train_atac_activity: ad.AnnData,
    val_query_atac: ad.AnnData,
    val_true_rna: ad.AnnData,
    val_atac_activity: ad.AnnData,
) -> None:
    safe_write_h5ad(train_rna_ref, outdir / "train_rna_ref.h5ad")
    safe_write_h5ad(train_atac_full, outdir / "train_atac_full.h5ad")
    safe_write_h5ad(train_atac_activity, outdir / "train_atac_activity.h5ad")
    safe_write_h5ad(val_query_atac, outdir / "val_query_atac.h5ad")
    safe_write_h5ad(val_true_rna, outdir / "val_true_rna.h5ad")
    safe_write_h5ad(val_atac_activity, outdir / "val_atac_activity.h5ad")


def generate_splits(
    rna_qc_path: Path,
    atac_qc_path: Path,
    atac_gas_path: Path,
    out_root: Path,
    cfg: Dict,
) -> pd.DataFrame:
    split_cfg = cfg["split"]
    checks = cfg["checks"]

    rna = sc.read_h5ad(str(rna_qc_path))
    atac_peak = sc.read_h5ad(str(atac_qc_path))
    atac_gas = sc.read_h5ad(str(atac_gas_path))

    common_cells = get_common_names(
        rna.obs_names.tolist(),
        atac_peak.obs_names.tolist(),
        atac_gas.obs_names.tolist(),
    )
    if len(common_cells) == 0:
        raise ValueError(
            "No common cells shared by RNA_counts_qc, "
            "ATAC_counts_qc, and ATAC_gas."
        )

    rna = rna[common_cells].copy()
    atac_peak = atac_peak[common_cells].copy()
    atac_gas = atac_gas[common_cells].copy()

    common_features = get_common_names(
        rna.var_names.tolist(),
        atac_gas.var_names.tolist(),
    )
    min_common = int(checks["min_common_features"])
    if len(common_features) < min_common:
        raise ValueError(
            "Too few common gene-level features between RNA and ATAC_gas: "
            f"{len(common_features)} < {min_common}"
        )

    # Retain the notebook behavior: splits use all common RNA/gene-activity
    # features, not only feature_aligned HVGs.
    rna = rna[:, common_features].copy()
    atac_gas = atac_gas[:, common_features].copy()

    all_cells = np.array(common_cells, dtype=object)
    seed = int(split_cfg["seed"])
    train_cells, val_cells = train_test_split(
        all_cells,
        test_size=float(split_cfg["val_frac"]),
        random_state=seed,
        shuffle=True,
    )
    train_cells = np.array(sorted(train_cells.tolist()), dtype=object)
    val_cells = np.array(sorted(val_cells.tolist()), dtype=object)

    ensure_dir(out_root)
    summaries = []

    for sf in [float(x) for x in split_cfg["single_fracs"]]:
        ratio_label = f"single_{int(round(sf * 100)):03d}"
        outdir = ensure_dir(out_root / ratio_label)

        train_assign = assign_partial_modality(
            train_cells,
            single_frac=sf,
            rna_keep_prob=float(split_cfg["rna_keep_prob"]),
            seed=seed + int(round(sf * 1000)) + 11,
        )
        val_assign = assign_partial_modality(
            val_cells,
            single_frac=sf,
            rna_keep_prob=float(split_cfg["rna_keep_prob"]),
            seed=seed + int(round(sf * 1000)) + 97,
        )

        # Retain original R-style logic.
        train_rna_ref_cells = sorted(
            train_assign.paired_cells + train_assign.rna_only_cells
        )
        train_atac_ref_cells = sorted(
            train_assign.paired_cells + train_assign.atac_only_cells
        )
        val_query_atac_cells = sorted(
            val_assign.paired_cells + val_assign.atac_only_cells
        )

        if (
            len(train_rna_ref_cells) < int(checks["min_train_rna_cells"])
            or len(train_atac_ref_cells) < int(checks["min_train_atac_cells"])
            or len(val_query_atac_cells) < int(checks["min_val_query_cells"])
        ):
            meta = {
                "ratio_label": ratio_label,
                "single_frac": sf,
                "status": "skipped",
                "reason": "too_few_cells",
                "train_rna_ref": len(train_rna_ref_cells),
                "train_atac_ref": len(train_atac_ref_cells),
                "val_query_atac": len(val_query_atac_cells),
            }
            with (outdir / "split_info.json").open(
                "w", encoding="utf-8"
            ) as f:
                json.dump(meta, f, indent=2, ensure_ascii=False)
            summaries.append(meta)
            continue

        train_rna_ref = subset_and_copy(rna, train_rna_ref_cells)

        # Retain the original naming/behavior:
        # train_atac_full contains ALL training cells at peak level.
        train_atac_full = subset_and_copy(atac_peak, train_cells)

        # Gene-activity reference uses paired + ATAC-only training cells.
        train_atac_activity = subset_and_copy(
            atac_gas,
            train_atac_ref_cells,
        )

        # Validation query is peak-level ATAC.
        val_query_atac = subset_and_copy(
            atac_peak,
            val_query_atac_cells,
        )

        # Ground-truth RNA and gene-activity are saved for the same query cells.
        val_true_rna = subset_and_copy(
            rna,
            val_query_atac_cells,
        )
        val_atac_activity = subset_and_copy(
            atac_gas,
            val_query_atac_cells,
        )

        save_split_h5ads(
            outdir,
            train_rna_ref=train_rna_ref,
            train_atac_full=train_atac_full,
            train_atac_activity=train_atac_activity,
            val_query_atac=val_query_atac,
            val_true_rna=val_true_rna,
            val_atac_activity=val_atac_activity,
        )

        meta = {
            "ratio_label": ratio_label,
            "single_frac": sf,
            "status": "ok",
            "train_cells": train_cells.tolist(),
            "val_cells": val_cells.tolist(),
            "train_paired_cells": train_assign.paired_cells,
            "train_rna_only_cells": train_assign.rna_only_cells,
            "train_atac_only_cells": train_assign.atac_only_cells,
            "val_paired_cells": val_assign.paired_cells,
            "val_rna_only_cells": val_assign.rna_only_cells,
            "val_atac_only_cells": val_assign.atac_only_cells,
            "train_rna_ref_cells": train_rna_ref_cells,
            "train_atac_ref_cells": train_atac_ref_cells,
            "val_query_atac_cells": val_query_atac_cells,
            "common_gene_features": common_features,
            "n_train_rna_ref_cells": len(train_rna_ref_cells),
            "n_train_atac_ref_cells": len(train_atac_ref_cells),
            "n_val_query_atac_cells": len(val_query_atac_cells),
            "n_common_features": len(common_features),
        }
        with (outdir / "split_info.json").open(
            "w", encoding="utf-8"
        ) as f:
            json.dump(meta, f, indent=2, ensure_ascii=False)

        summaries.append({
            "ratio_label": ratio_label,
            "single_frac": sf,
            "status": "ok",
            "train_paired": len(train_assign.paired_cells),
            "train_rna_only": len(train_assign.rna_only_cells),
            "train_atac_only": len(train_assign.atac_only_cells),
            "val_paired": len(val_assign.paired_cells),
            "val_rna_only": len(val_assign.rna_only_cells),
            "val_atac_only": len(val_assign.atac_only_cells),
            "train_rna_ref": len(train_rna_ref_cells),
            "train_atac_ref": len(train_atac_ref_cells),
            "val_query_atac": len(val_query_atac_cells),
            "n_common_features": len(common_features),
        })

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(out_root / "summary_all_ratios.csv", index=False)
    return summary_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="ATLAS PBMC RNA-ATAC preprocessing"
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to preprocessing/configs/pbmc.yaml",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    project_root = get_project_root(cfg)

    paths_cfg = cfg["paths"]
    input_h5 = require_file(
        resolve_project_path(project_root, paths_cfg["input_h5"]),
        "PBMC 10x Multiome H5",
    )
    gtf_path = require_file(
        resolve_project_path(project_root, paths_cfg["gtf"]),
        "GENCODE GTF",
    )

    fragments_path = resolve_project_path(
        project_root,
        paths_cfg.get("fragments"),
    )
    if fragments_path is not None and not fragments_path.exists():
        print(
            "Note: ATAC fragments file is not present. "
            "This is OK for this preprocessing pipeline because the original "
            "notebook only kept that path for compatibility."
        )

    interim_dir = ensure_dir(
        resolve_project_path(project_root, paths_cfg["interim"])
    )
    processed_dir = ensure_dir(
        resolve_project_path(project_root, paths_cfg["processed"])
    )

    rna_qc_path = interim_dir / "RNA_counts_qc.h5ad"
    atac_qc_path = interim_dir / "ATAC_counts_qc.h5ad"
    atac_gas_path = interim_dir / "ATAC_gas.h5ad"
    feature_aligned_path = interim_dir / "feature_aligned.h5ad"

    set_seed(int(cfg["split"]["seed"]))

    print("=" * 80)
    print("ATLAS PBMC RNA-ATAC preprocessing")
    print("=" * 80)
    print("Input H5:", input_h5)
    print("GTF:", gtf_path)
    print("Interim:", interim_dir)
    print("Processed:", processed_dir)
    if fragments_path is not None:
        print("Fragments (not consumed by this pipeline):", fragments_path)

    print("\n[1/6] Reading 10x Multiome H5...")
    rna_raw, atac_raw = read_10x_multiome_h5(input_h5)
    print("RNA raw shape:", rna_raw.shape)
    print("ATAC raw shape:", atac_raw.shape)

    print("\n[2/6] Preprocessing RNA...")
    rna_qc = preprocess_rna(rna_raw, cfg)
    safe_write_h5ad(rna_qc, rna_qc_path)
    print("Saved:", rna_qc_path, rna_qc.shape)

    print("\n[3/6] Preprocessing ATAC...")
    atac_qc = preprocess_atac(atac_raw, cfg)
    safe_write_h5ad(atac_qc, atac_qc_path)
    print("Saved:", atac_qc_path, atac_qc.shape)

    print("\n[4/6] Building ATAC gene activity...")
    atac_gas = build_gene_activity(
        atac_qc_path=atac_qc_path,
        gtf_path=gtf_path,
        cfg=cfg,
    )
    safe_write_h5ad(atac_gas, atac_gas_path)
    print("Saved:", atac_gas_path, atac_gas.shape)

    print("\n[5/6] Building feature-aligned RNA/ATAC-activity object...")
    feature_aligned = build_feature_aligned(
        rna_qc_path=rna_qc_path,
        atac_gas_path=atac_gas_path,
    )
    safe_write_h5ad(feature_aligned, feature_aligned_path)
    print("Saved:", feature_aligned_path, feature_aligned.shape)

    print("\n[6/6] Generating ratio split datasets...")
    split_root = ensure_dir(
        processed_dir / str(paths_cfg["split_dir_name"])
    )
    summary = generate_splits(
        rna_qc_path=rna_qc_path,
        atac_qc_path=atac_qc_path,
        atac_gas_path=atac_gas_path,
        out_root=split_root,
        cfg=cfg,
    )
    print("\nSummary:")
    print(summary)
    print("\nSaved:", split_root / "summary_all_ratios.csv")
    print("\nPBMC preprocessing finished successfully.")


if __name__ == "__main__":
    main()
