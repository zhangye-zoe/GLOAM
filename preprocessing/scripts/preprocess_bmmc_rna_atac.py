#!/usr/bin/env python3
"""
ATLAS BMMC RNA-ATAC preprocessing pipeline.

The complete original workflow is retained:
1. read the raw BMMC multiome h5ad;
2. split GEX and ATAC modalities;
3. preprocess/QC RNA and save RNA_counts_qc.h5ad;
4. preprocess/QC ATAC and save ATAC_counts_qc.h5ad;
5. build ATAC gene activity from the GTF and save ATAC_gas.h5ad;
6. build feature_aligned_rna_atac.h5ad;
7. generate train/validation datasets at multiple unpaired-data ratios.

The three core intermediate files are ALWAYS regenerated before split generation.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import anndata as ad
import episcanpy as epi
import numpy as np
import pandas as pd
import scanpy as sc
from muon import atac as ac
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


def read_bmmc_multiome(input_h5ad: Path) -> Tuple[ad.AnnData, ad.AnnData]:
    adata = sc.read_h5ad(str(input_h5ad))
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()

    if "feature_types" not in adata.var.columns:
        raise ValueError("feature_types not found in adata.var")

    feature_types = adata.var["feature_types"].astype(str)
    rna = adata[:, feature_types.eq("GEX").values].copy()
    atac = adata[:, feature_types.eq("ATAC").values].copy()

    if rna.n_vars == 0 or atac.n_vars == 0:
        raise ValueError("feature_types exists, but GEX or ATAC features were not found.")

    # Keep only modality-relevant embeddings, matching the original notebook.
    rna.obsm = {
        k: rna.obsm[k]
        for k in ["GEX_X_pca", "GEX_X_umap"]
        if k in rna.obsm
    }
    atac.obsm = {
        k: atac.obsm[k]
        for k in ["ATAC_gene_activity", "ATAC_lsi_full", "ATAC_lsi_red", "ATAC_umap"]
        if k in atac.obsm
    }
    return rna, atac


def preprocess_rna(rna: ad.AnnData, cfg: Dict) -> ad.AnnData:
    p = cfg["rna"]
    rna = rna.copy()
    batch_key = p["batch_key"]
    require_obs_column(rna, batch_key, "batch0")

    rna.var["mt"] = rna.var_names.astype(str).str.startswith(p["mt_prefix"])
    sc.pp.calculate_qc_metrics(
        rna,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    sc.pp.filter_genes(rna, min_cells=int(p["min_cells_per_gene"]))
    sc.pp.filter_cells(rna, min_genes=int(p["min_genes"]))
    rna = rna[rna.obs["n_genes_by_counts"] <= int(p["max_genes"]), :].copy()
    rna = rna[rna.obs["total_counts"] <= int(p["max_counts"]), :].copy()
    rna = rna[rna.obs["pct_counts_mt"] < float(p["max_mt"]), :].copy()

    rna.layers["counts"] = rna.X.copy()
    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(
        rna,
        batch_key=batch_key,
        min_mean=float(p["hvg_min_mean"]),
        max_mean=float(p["hvg_max_mean"]),
        min_disp=float(p["hvg_min_disp"]),
    )
    return rna


def preprocess_atac(atac: ad.AnnData, cfg: Dict) -> ad.AnnData:
    p = cfg["atac"]
    batch_key = cfg["rna"]["batch_key"]
    atac = atac.copy()
    require_obs_column(atac, batch_key, "batch0")

    sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_genes(atac, min_cells=int(p["min_cells_per_peak"]))
    sc.pp.filter_cells(atac, min_genes=int(p["min_features_per_cell"]))
    atac = atac[
        atac.obs["n_genes_by_counts"] <= int(p["max_features_per_cell"]), :
    ].copy()

    atac.layers["counts"] = atac.X.copy()
    ac.pp.tfidf(atac, scale_factor=float(p["tfidf_scale_factor"]))
    sc.pp.highly_variable_genes(
        atac,
        min_mean=float(p["hvg_min_mean"]),
        max_mean=float(p["hvg_max_mean"]),
        min_disp=float(p["hvg_min_disp"]),
        batch_key=batch_key,
    )

    # EpiScanpy geneactivity expects genomic intervals such as chr1:100-200.
    atac.var_names = [
        name.replace("-", ":", 1) if "-" in name and ":" not in name else name
        for name in atac.var_names.astype(str)
    ]
    return atac


def build_atac_gene_activity(atac_qc: ad.AnnData, gtf_path: Path, cfg: Dict) -> ad.AnnData:
    p = cfg["atac"]
    atac = atac_qc.copy()
    atac.X = atac.layers["counts"].copy()

    gas = epi.tl.geneactivity(
        atac,
        str(gtf_path),
        annotation=str(p["gene_activity_annotation"]),
    )
    if gas is None:
        # Some EpiScanpy versions modify the object in place.
        gas = atac

    gas = gas[:, ~gas.var_names.duplicated()].copy()
    gas.layers["counts"] = gas.X.copy()

    # Retain the original optional TF-IDF normalization behavior.
    try:
        ac.pp.tfidf(gas, scale_factor=float(p["tfidf_scale_factor"]))
    except Exception as exc:
        print(f"Warning: optional TF-IDF on ATAC gene activity was skipped: {exc!r}")

    sc.pp.highly_variable_genes(gas)
    return gas


def build_feature_aligned(rna_qc: ad.AnnData, atac_gas: ad.AnnData) -> ad.AnnData:
    rna_hvg = list(rna_qc.var_names[rna_qc.var["highly_variable"]].astype(str))
    gas_hvg = list(atac_gas.var_names[atac_gas.var["highly_variable"]].astype(str))
    union_genes = sorted(set(rna_hvg) | set(gas_hvg))

    rna = rna_qc.copy()
    gas = atac_gas.copy()
    combined = ad.concat(
        [rna, gas],
        join="outer",
        label="modality",
        keys=["rna", "atac_activity"],
    )
    combined.uns["rna_hvg"] = rna_hvg
    combined.uns["atac_hvg"] = gas_hvg
    combined.uns["rna_nz"] = sorted(set(union_genes) & set(rna.var_names.astype(str)))
    combined.uns["atac_nz"] = sorted(set(union_genes) & set(gas.var_names.astype(str)))
    combined = combined[:, union_genes].copy()
    return combined


def assign_partial_modality(
    cells: Sequence[str],
    single_frac: float = 0.2,
    rna_keep_prob: float = 0.5,
    seed: Optional[int] = None,
) -> Dict[str, object]:
    """Create paired / RNA-only / ATAC-only cell assignments."""
    rng = random.Random(seed) if seed is not None else random
    cells = list(cells)
    n_single = round(len(cells) * single_frac)
    single_cells = rng.sample(cells, n_single) if n_single > 0 else []
    paired_cells = sorted(set(cells) - set(single_cells))

    modality = {c: "paired" for c in cells}
    for cell in single_cells:
        modality[cell] = "RNA" if rng.random() < rna_keep_prob else "ATAC"

    return {
        "modality": modality,
        "paired_cells": paired_cells,
        "rna_only_cells": sorted([c for c, m in modality.items() if m == "RNA"]),
        "atac_only_cells": sorted([c for c, m in modality.items() if m == "ATAC"]),
    }


def save_split_h5ads(
    outdir: Path,
    train_rna_ref: ad.AnnData,
    train_atac_full: ad.AnnData,
    train_atac_activity: ad.AnnData,
    val_query_atac: ad.AnnData,
    val_true_rna: ad.AnnData,
    val_atac_activity: ad.AnnData,
) -> None:
    ensure_dir(outdir)
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
    check_cfg = cfg["checks"]

    rna = sc.read_h5ad(str(rna_qc_path))
    atac_peak = sc.read_h5ad(str(atac_qc_path))
    atac_gas = sc.read_h5ad(str(atac_gas_path))

    common_cells = get_common_names(
        rna.obs_names.tolist(),
        atac_peak.obs_names.tolist(),
        atac_gas.obs_names.tolist(),
    )
    if not common_cells:
        raise ValueError(
            "No common cells shared by RNA_counts_qc, ATAC_counts_qc, and ATAC_gas."
        )

    rna = rna[common_cells].copy()
    atac_peak = atac_peak[common_cells].copy()
    atac_gas = atac_gas[common_cells].copy()

    common_features = get_common_names(rna.var_names.tolist(), atac_gas.var_names.tolist())
    min_common_features = int(check_cfg["min_common_features"])
    if len(common_features) < min_common_features:
        raise ValueError(
            f"Too few common gene-level features between RNA and ATAC_gas: "
            f"{len(common_features)} < {min_common_features}"
        )

    rna = rna[:, common_features].copy()
    atac_gas = atac_gas[:, common_features].copy()

    seed = int(split_cfg["seed"])
    val_frac = float(split_cfg["val_frac"])
    single_fracs = [float(x) for x in split_cfg["single_fracs"]]
    rna_keep_prob = float(split_cfg["rna_keep_prob"])

    all_cells = np.array(common_cells, dtype=object)
    train_cells, val_cells = train_test_split(
        all_cells,
        test_size=val_frac,
        random_state=seed,
        shuffle=True,
    )
    train_cells = np.array(sorted(train_cells.tolist()), dtype=object)
    val_cells = np.array(sorted(val_cells.tolist()), dtype=object)

    ensure_dir(out_root)
    summaries = []

    for sf in single_fracs:
        ratio_label = f"single_{int(round(sf * 100)):03d}"
        outdir = ensure_dir(out_root / ratio_label)

        # Preserve the original notebook behavior: assignments consume the globally
        # seeded Python RNG sequentially across ratios.
        train_assign = assign_partial_modality(
            train_cells,
            single_frac=sf,
            rna_keep_prob=rna_keep_prob,
        )
        val_assign = assign_partial_modality(
            val_cells,
            single_frac=sf,
            rna_keep_prob=rna_keep_prob,
        )

        train_rna_ref_cells = sorted(
            train_assign["paired_cells"] + train_assign["rna_only_cells"]
        )
        train_atac_ref_cells = sorted(
            train_assign["paired_cells"] + train_assign["atac_only_cells"]
        )
        val_query_atac_cells = sorted(
            val_assign["paired_cells"] + val_assign["atac_only_cells"]
        )

        min_train_rna = int(check_cfg["min_train_rna_cells"])
        min_val_query = int(check_cfg["min_val_query_cells"])
        if (
            len(train_rna_ref_cells) < min_train_rna
            or len(train_atac_ref_cells) < min_train_rna
            or len(val_query_atac_cells) < min_val_query
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
            with (outdir / "split_info.json").open("w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2, ensure_ascii=False)
            summaries.append(meta)
            continue

        train_rna_ref = subset_and_copy(rna, train_rna_ref_cells)
        train_atac_full = subset_and_copy(atac_peak, train_cells)
        train_atac_activity = subset_and_copy(atac_gas, train_atac_ref_cells)
        val_query_atac = subset_and_copy(atac_peak, val_query_atac_cells)
        val_true_rna = subset_and_copy(rna, val_query_atac_cells)
        val_atac_activity = subset_and_copy(atac_gas, val_query_atac_cells)

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
            "train_paired_cells": train_assign["paired_cells"],
            "train_rna_only_cells": train_assign["rna_only_cells"],
            "train_atac_only_cells": train_assign["atac_only_cells"],
            "val_paired_cells": val_assign["paired_cells"],
            "val_rna_only_cells": val_assign["rna_only_cells"],
            "val_atac_only_cells": val_assign["atac_only_cells"],
            "train_rna_ref_cells": train_rna_ref_cells,
            "train_atac_ref_cells": train_atac_ref_cells,
            "val_query_atac_cells": val_query_atac_cells,
            "common_gene_features": common_features,
            "n_train_rna_ref_cells": len(train_rna_ref_cells),
            "n_train_atac_ref_cells": len(train_atac_ref_cells),
            "n_val_query_atac_cells": len(val_query_atac_cells),
            "n_common_features": len(common_features),
        }
        with (outdir / "split_info.json").open("w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2, ensure_ascii=False)

        summaries.append(
            {
                "ratio_label": ratio_label,
                "single_frac": sf,
                "status": "ok",
                "train_paired": len(train_assign["paired_cells"]),
                "train_rna_only": len(train_assign["rna_only_cells"]),
                "train_atac_only": len(train_assign["atac_only_cells"]),
                "val_paired": len(val_assign["paired_cells"]),
                "val_rna_only": len(val_assign["rna_only_cells"]),
                "val_atac_only": len(val_assign["atac_only_cells"]),
                "train_rna_ref": len(train_rna_ref_cells),
                "train_atac_ref": len(train_atac_ref_cells),
                "val_query_atac": len(val_query_atac_cells),
                "n_common_features": len(common_features),
            }
        )

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(out_root / "summary_all_ratios.csv", index=False)
    return summary_df


def run_pipeline(config_path: Path) -> pd.DataFrame:
    cfg = load_config(config_path)
    root = get_project_root(cfg)
    paths = cfg["paths"]

    input_h5ad = require_file(
        resolve_project_path(root, paths["multiome_input"]),
        "BMMC multiome input",
    )
    gtf_path = require_file(
        resolve_project_path(root, paths["gtf"]),
        "GENCODE GTF annotation",
    )
    interim_dir = ensure_dir(resolve_project_path(root, paths["rna_atac_interim"]))
    processed_dir = ensure_dir(resolve_project_path(root, paths["rna_atac_processed"]))

    seed = int(cfg["split"]["seed"])
    set_seed(seed)

    rna_qc_path = interim_dir / "RNA_counts_qc.h5ad"
    atac_qc_path = interim_dir / "ATAC_counts_qc.h5ad"
    atac_gas_path = interim_dir / "ATAC_gas.h5ad"
    feature_aligned_path = interim_dir / "feature_aligned_rna_atac.h5ad"
    split_root = processed_dir / "results_ratio_loop_rna_atac"

    print(f"Project root: {root}")
    print(f"Input: {input_h5ad}")
    print("Reading BMMC multiome...")
    rna_raw, atac_raw = read_bmmc_multiome(input_h5ad)
    print("RNA raw:", rna_raw.shape)
    print("ATAC raw:", atac_raw.shape)

    print("[1/5] Preprocess RNA -> RNA_counts_qc.h5ad")
    rna_qc = preprocess_rna(rna_raw, cfg)
    print("RNA QC:", rna_qc.shape)
    safe_write_h5ad(rna_qc, rna_qc_path)

    print("[2/5] Preprocess ATAC -> ATAC_counts_qc.h5ad")
    atac_qc = preprocess_atac(atac_raw, cfg)
    print("ATAC QC:", atac_qc.shape)
    safe_write_h5ad(atac_qc, atac_qc_path)

    print("[3/5] Build ATAC gene activity -> ATAC_gas.h5ad")
    atac_gas = build_atac_gene_activity(atac_qc, gtf_path, cfg)
    print("ATAC gene activity:", atac_gas.shape)
    safe_write_h5ad(atac_gas, atac_gas_path)

    print("[4/5] Build feature-aligned RNA-ATAC file")
    feature_aligned = build_feature_aligned(rna_qc, atac_gas)
    print("Feature aligned:", feature_aligned.shape)
    safe_write_h5ad(feature_aligned, feature_aligned_path)

    print("[5/5] Generate ratio splits")
    summary_df = generate_splits(
        rna_qc_path=rna_qc_path,
        atac_qc_path=atac_qc_path,
        atac_gas_path=atac_gas_path,
        out_root=split_root,
        cfg=cfg,
    )
    print(summary_df.to_string(index=False))
    print(f"RNA-ATAC preprocessing complete. Outputs: {interim_dir} and {processed_dir}")
    return summary_df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Preprocess BMMC RNA-ATAC data for ATLAS.")
    default_config = Path(__file__).resolve().parents[1] / "configs" / "bmmc.yaml"
    parser.add_argument(
        "--config",
        type=Path,
        default=default_config,
        help=f"YAML config path (default: {default_config})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_pipeline(args.config)


if __name__ == "__main__":
    main()
