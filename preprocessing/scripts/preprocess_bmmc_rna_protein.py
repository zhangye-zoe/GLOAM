#!/usr/bin/env python3
"""
ATLAS BMMC RNA-Protein / CITE-seq preprocessing pipeline.

Complete workflow retained from the original notebook:
1. read the raw BMMC CITE-seq h5ad;
2. split GEX and ADT modalities;
3. preprocess/QC RNA;
4. preprocess Protein/ADT with CLR normalization;
5. optionally map ADT marker names to HGNC-approved gene symbols;
6. add the original coarse BMMC cell-type grouping;
7. write intermediate h5ad files;
8. generate train/validation datasets at multiple unpaired-data ratios.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.model_selection import train_test_split

try:
    from muon import prot as pt
except Exception as exc:
    pt = None
    print("Warning: failed to import muon.prot. Protein CLR preprocessing will fail.")
    print(repr(exc))

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


# Original coarse BMMC label grouping retained from the notebook.
COARSE_CELLTYPE_MAPPING = {
    # CD4 T
    "CD4+ T naive": "CD4 T",
    "CD4+ T activated": "CD4 T",
    "CD4+ T activated integrinB7+": "CD4 T",
    "CD4+ T CD314+ CD45RA+": "CD4 T",
    "T reg": "CD4 T",

    # CD8 T
    "CD8+ T naive": "CD8 T",
    "CD8+ T naive CD127+ CD26- CD101-": "CD8 T",
    "CD8+ T CD49f+": "CD8 T",
    "CD8+ T TIGIT+ CD45RO+": "CD8 T",
    "CD8+ T CD57+ CD45RA+": "CD8 T",
    "CD8+ T CD69+ CD45RO+": "CD8 T",
    "CD8+ T TIGIT+ CD45RA+": "CD8 T",
    "CD8+ T CD69+ CD45RA+": "CD8 T",
    "CD8+ T CD57+ CD45RO+": "CD8 T",
    "CD8+ T": "CD8 T",
    "MAIT": "CD8 T",
    "gdT TCRVD2+": "CD8 T",
    "gdT CD158b+": "CD8 T",
    "dnT": "CD8 T",

    # B cells
    "Naive CD20+ B IGKC+": "B cell",
    "Naive CD20+ B IGKC-": "B cell",
    "Naive CD20+ B": "B cell",
    "B1 B IGKC+": "B cell",
    "B1 B IGKC-": "B cell",
    "B1 B": "B cell",
    "Transitional B": "B cell",

    # Plasma cells
    "Plasmablast IGKC+": "Plasma cell",
    "Plasmablast IGKC-": "Plasma cell",
    "Plasma cell IGKC+": "Plasma cell",
    "Plasma cell IGKC-": "Plasma cell",
    "Plasma cell": "Plasma cell",

    # NK
    "NK": "NK",
    "NK CD158e1+": "NK",

    # Monocytes
    "CD14+ Mono": "Mono",
    "CD16+ Mono": "Mono",

    # DC
    "pDC": "DC",
    "cDC1": "DC",
    "cDC2": "DC",

    # ILC
    "ILC": "ILC",
    "ILC1": "ILC",

    # Progenitors
    "HSC": "Progenitor",
    "Lymph prog": "Progenitor",
    "G/M prog": "Progenitor",
    "MK/E prog": "Progenitor",
    "ID2-hi myeloid prog": "Progenitor",
    "T prog cycling": "Progenitor",

    # Erythroid
    "Erythroblast": "Erythroid",
    "Normoblast": "Erythroid",
    "Proerythroblast": "Erythroid",
    "Reticulocyte": "Erythroid",
}


def set_seed(seed: int = 1234) -> None:
    random.seed(seed)
    np.random.seed(seed)


def add_coarse_celltype(
    adata: ad.AnnData,
    source_key: str,
    target_key: str = "celltype",
) -> None:
    """Add the original coarse BMMC label while preserving the original cell_type column."""
    if source_key in adata.obs.columns:
        adata.obs[target_key] = adata.obs[source_key].map(COARSE_CELLTYPE_MAPPING)
        adata.obs[target_key] = adata.obs[target_key].fillna(
            adata.obs[source_key].astype(str)
        )


def read_bmmc_cite(input_h5ad: Path) -> Tuple[ad.AnnData, ad.AnnData]:
    adata = sc.read_h5ad(str(input_h5ad))
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()

    if "cellid" not in adata.obs.columns:
        adata.obs["cellid"] = np.arange(adata.n_obs)

    if "feature_types" not in adata.var.columns:
        raise ValueError(
            "feature_types not found in adata.var. Please add feature_types or manually split RNA/Protein."
        )

    feature_types = adata.var["feature_types"].astype(str)
    gex_mask = feature_types.eq("GEX")
    adt_mask = feature_types.eq("ADT")
    if gex_mask.sum() == 0 or adt_mask.sum() == 0:
        raise ValueError("feature_types exists, but GEX or ADT features were not found.")

    rna = adata[:, gex_mask.values].copy()
    prot = adata[:, adt_mask.values].copy()

    rna.obsm = {
        k: rna.obsm[k]
        for k in ["GEX_X_pca", "GEX_X_umap"]
        if k in rna.obsm
    }
    prot.obsm = {
        k: prot.obsm[k]
        for k in ["ADT_X_pca", "ADT_X_umap", "ADT_isotype_controls"]
        if k in prot.obsm
    }

    rna.var_names_make_unique()
    prot.var_names_make_unique()
    return rna, prot


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

    add_coarse_celltype(rna, source_key=p["cell_type_key"])
    return rna


def load_hgnc_cd_mapper(mapping_csv: Optional[Path]):
    """Return a function mapping ADT marker names to HGNC-approved symbols."""
    if mapping_csv is None:
        return None
    if not mapping_csv.exists():
        print(f"HGNC mapping file not found; skip protein marker mapping: {mapping_csv}")
        return None

    df_hgnc = pd.read_csv(mapping_csv)
    required = {"Approved symbol", "Previous symbols", "Aliases"}
    missing = required - set(df_hgnc.columns)
    if missing:
        print(f"HGNC mapping file misses columns {missing}; skip protein marker mapping.")
        return None

    df_hgnc[["Previous symbols", "Aliases"]] = df_hgnc[
        ["Previous symbols", "Aliases"]
    ].fillna("")

    lookup = {}
    for _, row in df_hgnc.iterrows():
        approved = str(row["Approved symbol"]).strip()
        if not approved or approved == "nan":
            continue
        candidates = [approved]
        candidates += [
            s.strip()
            for s in str(row["Previous symbols"]).split(",")
            if s.strip()
        ]
        candidates += [
            s.strip()
            for s in str(row["Aliases"]).split(",")
            if s.strip()
        ]
        for candidate in candidates:
            lookup.setdefault(candidate, approved)

    def mapper(marker_name: str) -> str:
        return lookup.get(str(marker_name), str(marker_name))

    return mapper


def map_protein_var_names(
    prot: ad.AnnData,
    mapping_csv: Optional[Path],
) -> ad.AnnData:
    """Retain original ADT names and optionally standardize var_names to HGNC symbols."""
    prot = prot.copy()
    mapper = load_hgnc_cd_mapper(mapping_csv)

    prot.var["cd_name"] = prot.var_names.astype(str)
    if mapper is None:
        prot.var["gene_name"] = prot.var["cd_name"].astype(str)
    else:
        prot.var["gene_name"] = [mapper(x) for x in prot.var["cd_name"].astype(str)]

    prot.var_names = prot.var["gene_name"].astype(str).tolist()
    prot.var_names_make_unique()
    return prot


def preprocess_protein(
    prot: ad.AnnData,
    cfg: Dict,
    mapping_csv: Optional[Path],
) -> ad.AnnData:
    if pt is None:
        raise ImportError(
            "muon.prot is not available. Install muon before running Protein CLR preprocessing."
        )

    p = cfg["protein"]
    rna_cfg = cfg["rna"]
    prot = prot.copy()
    require_obs_column(prot, rna_cfg["batch_key"], "batch0")

    sc.pp.calculate_qc_metrics(prot, percent_top=None, log1p=False, inplace=True)
    min_cells = int(p["min_cells"])
    min_counts_cell = int(p["min_counts_cell"])
    if min_cells > 0:
        sc.pp.filter_genes(prot, min_cells=min_cells)
    if min_counts_cell > 0:
        prot = prot[prot.obs["total_counts"] >= min_counts_cell, :].copy()

    prot.layers["counts"] = prot.X.copy()

    # CLR normalization updates prot.X.
    pt.pp.clr(prot)

    # Preserve original downstream compatibility behavior.
    prot.var["highly_variable"] = True

    # Optional ADT/CD marker -> HGNC symbol mapping.
    prot = map_protein_var_names(prot, mapping_csv)

    add_coarse_celltype(prot, source_key=rna_cfg["cell_type_key"])
    return prot


def build_feature_aligned_rna_protein(
    rna_qc: ad.AnnData,
    prot_qc: ad.AnnData,
) -> ad.AnnData:
    """Build the original convenience combined RNA-Protein AnnData."""
    rna = rna_qc.copy()
    prot = prot_qc.copy()

    rna.var["modality_feature_type"] = "RNA"
    prot.var["modality_feature_type"] = "Protein"

    combined = ad.concat(
        [rna, prot],
        join="outer",
        label="modality",
        keys=["rna", "protein"],
        index_unique="__",
    )
    hv = (
        rna.var["highly_variable"]
        if "highly_variable" in rna.var.columns
        else np.ones(rna.n_vars, dtype=bool)
    )
    combined.uns["rna_hvg"] = list(rna.var_names[hv].astype(str))
    combined.uns["protein_features"] = list(prot.var_names.astype(str))
    return combined


def assign_partial_modality(
    cells: Sequence[str],
    single_frac: float = 0.2,
    rna_keep_prob: float = 0.5,
    seed: Optional[int] = None,
) -> Dict[str, object]:
    """Create paired / RNA-only / Protein-only cell assignments."""
    rng = random.Random(seed) if seed is not None else random
    cells = list(cells)
    n_single = round(len(cells) * single_frac)
    single_cells = rng.sample(cells, n_single) if n_single > 0 else []
    paired_cells = sorted(set(cells) - set(single_cells))

    modality = {c: "paired" for c in cells}
    for cell in single_cells:
        modality[cell] = "RNA" if rng.random() < rna_keep_prob else "Protein"

    return {
        "modality": modality,
        "paired_cells": paired_cells,
        "rna_only_cells": sorted([c for c, m in modality.items() if m == "RNA"]),
        "protein_only_cells": sorted(
            [c for c, m in modality.items() if m == "Protein"]
        ),
    }


def save_split_h5ads(
    outdir: Path,
    train_rna_ref: ad.AnnData,
    train_protein_ref: ad.AnnData,
    train_protein_full: ad.AnnData,
    val_query_protein: ad.AnnData,
    val_true_rna: ad.AnnData,
    val_true_protein: ad.AnnData,
) -> None:
    ensure_dir(outdir)
    safe_write_h5ad(train_rna_ref, outdir / "train_rna_ref.h5ad")
    safe_write_h5ad(train_protein_ref, outdir / "train_protein_ref.h5ad")
    safe_write_h5ad(train_protein_full, outdir / "train_protein_full.h5ad")
    safe_write_h5ad(val_query_protein, outdir / "val_query_protein.h5ad")
    safe_write_h5ad(val_true_rna, outdir / "val_true_rna.h5ad")
    safe_write_h5ad(val_true_protein, outdir / "val_true_protein.h5ad")


def generate_splits(
    rna_qc_path: Path,
    protein_qc_path: Path,
    out_root: Path,
    cfg: Dict,
) -> pd.DataFrame:
    split_cfg = cfg["split"]
    check_cfg = cfg["checks"]

    rna = sc.read_h5ad(str(rna_qc_path))
    prot = sc.read_h5ad(str(protein_qc_path))

    common_cells = get_common_names(rna.obs_names.tolist(), prot.obs_names.tolist())
    if not common_cells:
        raise ValueError("No common cells shared by RNA_counts_qc and protein_counts_qc.")

    rna = rna[common_cells].copy()
    prot = prot[common_cells].copy()

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

        train_assign = assign_partial_modality(
            train_cells,
            single_frac=sf,
            rna_keep_prob=rna_keep_prob,
            seed=seed + int(round(sf * 1000)) + 11,
        )
        val_assign = assign_partial_modality(
            val_cells,
            single_frac=sf,
            rna_keep_prob=rna_keep_prob,
            seed=seed + int(round(sf * 1000)) + 97,
        )

        train_rna_ref_cells = sorted(
            train_assign["paired_cells"] + train_assign["rna_only_cells"]
        )
        train_protein_ref_cells = sorted(
            train_assign["paired_cells"] + train_assign["protein_only_cells"]
        )
        val_query_protein_cells = sorted(
            val_assign["paired_cells"] + val_assign["protein_only_cells"]
        )

        if (
            len(train_rna_ref_cells) < int(check_cfg["min_train_rna_cells"])
            or len(train_protein_ref_cells) < int(check_cfg["min_train_protein_cells"])
            or len(val_query_protein_cells) < int(check_cfg["min_val_query_cells"])
        ):
            meta = {
                "ratio_label": ratio_label,
                "single_frac": sf,
                "status": "skipped",
                "reason": "too_few_cells",
                "train_rna_ref": len(train_rna_ref_cells),
                "train_protein_ref": len(train_protein_ref_cells),
                "val_query_protein": len(val_query_protein_cells),
            }
            with (outdir / "split_info.json").open("w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2, ensure_ascii=False)
            summaries.append(meta)
            continue

        train_rna_ref = subset_and_copy(rna, train_rna_ref_cells)
        train_protein_ref = subset_and_copy(prot, train_protein_ref_cells)
        train_protein_full = subset_and_copy(prot, train_cells)
        val_query_protein = subset_and_copy(prot, val_query_protein_cells)
        val_true_rna = subset_and_copy(rna, val_query_protein_cells)
        val_true_protein = subset_and_copy(prot, val_query_protein_cells)

        save_split_h5ads(
            outdir,
            train_rna_ref=train_rna_ref,
            train_protein_ref=train_protein_ref,
            train_protein_full=train_protein_full,
            val_query_protein=val_query_protein,
            val_true_rna=val_true_rna,
            val_true_protein=val_true_protein,
        )

        meta = {
            "ratio_label": ratio_label,
            "single_frac": sf,
            "status": "ok",
            "train_cells": train_cells.tolist(),
            "val_cells": val_cells.tolist(),
            "train_paired_cells": train_assign["paired_cells"],
            "train_rna_only_cells": train_assign["rna_only_cells"],
            "train_protein_only_cells": train_assign["protein_only_cells"],
            "val_paired_cells": val_assign["paired_cells"],
            "val_rna_only_cells": val_assign["rna_only_cells"],
            "val_protein_only_cells": val_assign["protein_only_cells"],
            "train_rna_ref_cells": train_rna_ref_cells,
            "train_protein_ref_cells": train_protein_ref_cells,
            "val_query_protein_cells": val_query_protein_cells,
            "rna_features": rna.var_names.astype(str).tolist(),
            "protein_features": prot.var_names.astype(str).tolist(),
            "n_train_rna_ref_cells": len(train_rna_ref_cells),
            "n_train_protein_ref_cells": len(train_protein_ref_cells),
            "n_val_query_protein_cells": len(val_query_protein_cells),
            "n_rna_features": rna.n_vars,
            "n_protein_features": prot.n_vars,
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
                "train_protein_only": len(train_assign["protein_only_cells"]),
                "val_paired": len(val_assign["paired_cells"]),
                "val_rna_only": len(val_assign["rna_only_cells"]),
                "val_protein_only": len(val_assign["protein_only_cells"]),
                "train_rna_ref": len(train_rna_ref_cells),
                "train_protein_ref": len(train_protein_ref_cells),
                "val_query_protein": len(val_query_protein_cells),
                "n_rna_features": rna.n_vars,
                "n_protein_features": prot.n_vars,
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
        resolve_project_path(root, paths["cite_input"]),
        "BMMC CITE-seq input",
    )

    mapping_value = paths.get("hgnc_cd_mapping")
    mapping_csv = resolve_project_path(root, mapping_value) if mapping_value else None
    if mapping_csv is not None and not mapping_csv.exists():
        print(
            "Warning: HGNC/CD mapping file is absent. The pipeline will continue without "
            f"protein marker renaming: {mapping_csv}"
        )

    interim_dir = ensure_dir(resolve_project_path(root, paths["rna_protein_interim"]))
    processed_dir = ensure_dir(resolve_project_path(root, paths["rna_protein_processed"]))

    seed = int(cfg["split"]["seed"])
    set_seed(seed)

    rna_qc_path = interim_dir / "RNA_counts_qc.h5ad"
    protein_qc_path = interim_dir / "protein_counts_qc.h5ad"
    feature_aligned_path = interim_dir / "feature_aligned_rna_protein.h5ad"
    split_root = processed_dir / "results_ratio_loop_rna_protein"

    print(f"Project root: {root}")
    print(f"Input: {input_h5ad}")
    print("Reading BMMC CITE-seq RNA/Protein data...")
    rna_raw, prot_raw = read_bmmc_cite(input_h5ad)
    print("RNA raw:", rna_raw.shape)
    print("Protein raw:", prot_raw.shape)

    print("[1/4] Preprocess RNA -> RNA_counts_qc.h5ad")
    rna_qc = preprocess_rna(rna_raw, cfg)
    print("RNA QC:", rna_qc.shape)
    safe_write_h5ad(rna_qc, rna_qc_path)

    print("[2/4] Preprocess Protein -> protein_counts_qc.h5ad")
    prot_qc = preprocess_protein(prot_raw, cfg, mapping_csv=mapping_csv)
    print("Protein QC:", prot_qc.shape)
    safe_write_h5ad(prot_qc, protein_qc_path)

    print("[3/4] Build feature-aligned RNA-Protein file")
    feature_aligned = build_feature_aligned_rna_protein(rna_qc, prot_qc)
    print("Feature aligned:", feature_aligned.shape)
    safe_write_h5ad(feature_aligned, feature_aligned_path)

    print("[4/4] Generate ratio splits")
    summary_df = generate_splits(
        rna_qc_path=rna_qc_path,
        protein_qc_path=protein_qc_path,
        out_root=split_root,
        cfg=cfg,
    )
    print(summary_df.to_string(index=False))
    print(f"RNA-Protein preprocessing complete. Outputs: {interim_dir} and {processed_dir}")
    return summary_df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Preprocess BMMC RNA-Protein/CITE-seq data for ATLAS."
    )
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
