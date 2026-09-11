#!/usr/bin/env python3
"""Shared utilities for ATLAS preprocessing scripts."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import anndata as ad
import numpy as np
import yaml


def load_config(config_path: os.PathLike) -> Dict[str, Any]:
    """Load a YAML configuration file."""
    config_path = Path(config_path).expanduser().resolve()

    if not config_path.exists():
        raise FileNotFoundError(
            f"Config file not found: {config_path}"
        )

    with config_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    if not isinstance(cfg, dict):
        raise ValueError(
            f"Invalid YAML config: {config_path}"
        )

    return cfg


def get_project_root(cfg: Dict[str, Any]) -> Path:
    """
    Return the project root directory.

    If `project_root` is specified in the YAML config, use it.
    Otherwise infer the ATLAS root from this file location:

        ATLAS/
        └── preprocessing/
            └── scripts/
                └── common.py
    """
    root = cfg.get("project_root")

    if root:
        return Path(root).expanduser().resolve()

    return Path(__file__).resolve().parents[2]


def resolve_project_path(
    project_root: Path,
    value: Optional[str],
) -> Optional[Path]:
    """
    Resolve a path relative to the project root.

    Absolute paths are returned unchanged.
    """
    if value is None:
        return None

    path = Path(value).expanduser()

    if path.is_absolute():
        return path

    return (project_root / path).resolve()


def ensure_dir(path: os.PathLike) -> Path:
    """Create a directory if it does not already exist."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)

    return path


def require_file(
    path: os.PathLike,
    label: str,
) -> Path:
    """Check that a required file exists."""
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(
            f"{label} not found: {path}"
        )

    return path


def safe_write_h5ad(
    adata: ad.AnnData,
    path: os.PathLike,
) -> None:
    """
    Safely write an AnnData object.

    Converts np.matrix objects to ndarray before serialization,
    because np.matrix objects may cause issues when writing H5AD files.
    """
    path = Path(path)
    ensure_dir(path.parent)

    obj = adata.copy()

    if isinstance(obj.X, np.matrix):
        obj.X = np.asarray(obj.X)

    for key in list(obj.layers.keys()):
        if isinstance(obj.layers[key], np.matrix):
            obj.layers[key] = np.asarray(
                obj.layers[key]
            )

    obj.write(str(path))


def require_obs_column(
    adata: ad.AnnData,
    key: str,
    fill_value: str = "batch0",
) -> None:
    """
    Ensure an AnnData.obs column exists.

    If it does not exist, create it using `fill_value`.
    """
    if key not in adata.obs.columns:
        adata.obs[key] = fill_value


def get_common_names(
    *arrays: Sequence[str],
) -> List[str]:
    """
    Return sorted names shared by all input sequences.
    """
    if not arrays:
        return []

    common = set(arrays[0])

    for arr in arrays[1:]:
        common &= set(arr)

    return sorted(common)


def subset_and_copy(
    adata: ad.AnnData,
    cells: Sequence[str],
    features: Optional[Sequence[str]] = None,
) -> ad.AnnData:
    """
    Subset an AnnData object by cells and, optionally, features.

    Always returns a copy.
    """
    out = adata[list(cells)].copy()

    if features is not None:
        out = out[:, list(features)].copy()

    return out