# ATLAS Data Preprocessing

This directory contains reproducible preprocessing pipelines for the **PBMC** and **BMMC** multi-omics datasets used in ATLAS.

The original preprocessing code was developed in Jupyter notebooks and has been reorganized into standalone Python scripts with YAML configuration files. The current structure separates:

- raw public datasets;
- external annotation resources;
- intermediate preprocessing outputs;
- final model-ready split datasets;
- preprocessing code, configuration files, and logs.

The pipelines currently cover:

- **PBMC RNA–ATAC**
- **BMMC RNA–ATAC**
- **BMMC RNA–Protein**

---

## 1. Directory layout

The expected project structure is:

```text
/data5/zhangye/ATLAS/
├── data/
│   ├── raw/
│   │   ├── PBMC/
│   │   │   ├── pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
│   │   │   └── pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
│   │   │
│   │   └── BMMC/
│   │       ├── GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad
│   │       └── GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad
│   │
│   ├── external/
│   │   └── annotations/
│   │       ├── gencode.v38.primary_assembly.annotation.gtf
│   │       └── hgnc_cd_gene_mapping.csv
│   │
│   ├── interim/
│   │   ├── PBMC/
│   │   │   └── RNA_ATAC/
│   │   └── BMMC/
│   │       ├── RNA_ATAC/
│   │       └── RNA_PROTEIN/
│   │
│   └── processed/
│       ├── PBMC/
│       │   └── RNA_ATAC/
│       └── BMMC/
│           ├── RNA_ATAC/
│           └── RNA_PROTEIN/
│
└── preprocessing/
    ├── README.md
    ├── requirements.txt
    ├── requirements_pbmc.txt
    ├── run_bmmc.sh
    ├── run_pbmc.sh
    ├── configs/
    │   ├── bmmc.yaml
    │   └── pbmc.yaml
    ├── logs/
    │   └── .gitkeep
    └── scripts/
        ├── common.py
        ├── preprocess_bmmc_rna_atac.py
        ├── preprocess_bmmc_rna_protein.py
        └── preprocess_pbmc_rna_atac.py
```

`common.py` contains shared utility functions used by both PBMC and BMMC preprocessing scripts.

---

## 2. Download the public datasets

### 2.1 PBMC RNA–ATAC dataset

Source:

https://www.10xgenomics.com/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-x-1-standard-2-0-0

This is a 10x Genomics **Single Cell Multiome ATAC + Gene Expression** dataset.

On the 10x Genomics dataset page, open **Output and supplemental files**.

For the current preprocessing pipeline, download:

#### Required

**Filtered feature barcode matrix (HDF5)**

This file contains both:

- Gene Expression features;
- ATAC Peak features.

The preprocessing script reads this HDF5 file with `scanpy.read_10x_h5(..., gex_only=False)` and separates the two modalities using `feature_types`.

The current `pbmc.yaml` expects the local filename:

```text
pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
```

Place it at:

```text
/data5/zhangye/ATLAS/data/raw/PBMC/
```

#### Optional for the current preprocessing pipeline

**ATAC Per fragment information file (TSV.GZ)**

The current local filename is:

```text
pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
```

Place it at:

```text
/data5/zhangye/ATLAS/data/raw/PBMC/
```

The fragments file is retained for provenance and future ATAC analyses. The current preprocessing script does **not** use it to generate the RNA, peak-level ATAC, or gene-activity matrices.

> **Note**
>
> 10x Genomics may use a different sample prefix for files downloaded from the dataset page. If the downloaded files have different names, either:
>
> 1. rename them to the filenames above; or
> 2. update the corresponding paths in `preprocessing/configs/pbmc.yaml`.

The file `pbmc_azimuth_annotations.csv` is not required by the current PBMC preprocessing pipeline.

---

### 2.2 BMMC dataset

Source:

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122

The BMMC dataset contains multi-omics measurements from bone marrow mononuclear cells and was released for the NeurIPS 2021 multimodal single-cell challenge.

Download these two processed H5AD files from the **Supplementary files** section of the GSE194122 page:

```text
GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz
GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz
```

They correspond to:

```text
multiome_BMMC
    RNA + ATAC

cite_BMMC
    RNA + Protein / ADT
```

After downloading, decompress them:

```bash
gunzip GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz
gunzip GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz
```

Then place the resulting files at:

```text
/data5/zhangye/ATLAS/data/raw/BMMC/
```

The final directory should contain:

```text
data/raw/BMMC/
├── GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad
└── GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad
```

The preprocessing scripts split the modalities from these H5AD objects and regenerate the downstream intermediate files.

---

## 3. External annotation files

The public datasets above are not the only files required by the preprocessing code.

External annotation resources are stored separately under:

```text
/data5/zhangye/ATLAS/data/external/annotations/
```

### 3.1 GENCODE gene annotation

Required file:

```text
gencode.v38.primary_assembly.annotation.gtf
```

Expected location:

```text
/data5/zhangye/ATLAS/data/external/annotations/gencode.v38.primary_assembly.annotation.gtf
```

This GTF is required by both RNA–ATAC pipelines:

```text
PBMC RNA–ATAC
BMMC RNA–ATAC
```

It is used by EpiScanpy to convert peak-level ATAC measurements into a gene-level **gene activity matrix**.

If the file already exists in an older project directory, copy it into ATLAS:

```bash
mkdir -p /data5/zhangye/ATLAS/data/external/annotations

cp /path/to/gencode.v38.primary_assembly.annotation.gtf \
   /data5/zhangye/ATLAS/data/external/annotations/
```

### 3.2 HGNC / CD marker mapping

Optional file:

```text
hgnc_cd_gene_mapping.csv
```

Expected location:

```text
/data5/zhangye/ATLAS/data/external/annotations/hgnc_cd_gene_mapping.csv
```

This file is used only by:

```text
BMMC RNA–Protein
```

It maps ADT/CD marker names to HGNC-approved gene symbols when a mapping is available.

For example, the mapping step can standardize protein feature names such as CD-marker aliases to gene symbols while retaining the original marker information.

The RNA–Protein pipeline can still run without this file; in that case, the original protein feature names are retained.

If using the original file from the previous project:

```bash
cp /data5/zhangye/scMRDR/scripts/BMMC/CD_gene/group-471.csv \
   /data5/zhangye/ATLAS/data/external/annotations/hgnc_cd_gene_mapping.csv
```

---

## 4. Preprocessing scripts

### `scripts/common.py`

Shared utilities used by all preprocessing pipelines, including:

- YAML configuration loading;
- project-relative path resolution;
- directory creation;
- required-file validation;
- safe AnnData serialization;
- shared cell/feature intersection;
- AnnData subsetting.

### `scripts/preprocess_pbmc_rna_atac.py`

Processes the 10x PBMC Multiome dataset.

Main workflow:

```text
10x filtered feature-barcode H5
        │
        ├── Gene Expression
        │       │
        │       ├── QC metrics
        │       ├── total-count normalization
        │       ├── log1p
        │       └── highly variable gene selection
        │
        └── ATAC Peaks
                │
                ├── QC metrics
                ├── TF-IDF
                └── highly variable peak selection
                        │
                        └── GENCODE GTF
                                │
                                v
                        ATAC gene activity
                                │
                                v
                        feature alignment
                                │
                                v
                    partial-pairing ratio splits
```

The PBMC preprocessing intentionally follows the original notebook behavior: QC metrics are calculated, but the explicit BMMC-style cell/feature filtering thresholds are not applied.

### `scripts/preprocess_bmmc_rna_atac.py`

Processes the BMMC Multiome H5AD file.

Main workflow:

```text
BMMC Multiome H5AD
        │
        ├── RNA preprocessing
        │
        └── ATAC preprocessing
                 │
                 └── GENCODE GTF
                         │
                         v
                 ATAC gene activity
                         │
                         v
                 feature alignment
                         │
                         v
             partial-pairing ratio splits
```

### `scripts/preprocess_bmmc_rna_protein.py`

Processes the BMMC CITE-seq H5AD file.

Main workflow:

```text
BMMC CITE-seq H5AD
        │
        ├── RNA preprocessing
        │
        └── Protein / ADT preprocessing
                 │
                 ├── CLR normalization
                 ├── optional HGNC/CD mapping
                 └── coarse cell-type annotation
                         │
                         v
                 partial-pairing ratio splits
```

The original fine-grained `obs["cell_type"]` annotation is preserved. A coarse annotation is added as:

```text
obs["celltype"]
```

---

## 5. Run preprocessing

Run all commands from the ATLAS project root:

```bash
cd /data5/zhangye/ATLAS
```

### 5.1 Install dependencies

If the existing ATLAS/scMRDR conda environment already contains the required packages, it is preferable to reuse that environment rather than upgrading packages unnecessarily.

For BMMC:

```bash
pip install -r preprocessing/requirements.txt
```

For PBMC-specific dependencies:

```bash
pip install -r preprocessing/requirements_pbmc.txt
```

Important packages include:

```text
anndata
scanpy
numpy
pandas
scipy
scikit-learn
PyYAML
h5py
episcanpy
muon
```

### 5.2 Run PBMC RNA–ATAC

Directly:

```bash
python preprocessing/scripts/preprocess_pbmc_rna_atac.py \
  --config preprocessing/configs/pbmc.yaml \
  2>&1 | tee preprocessing/logs/pbmc_rna_atac.log
```

Or use the wrapper:

```bash
bash preprocessing/run_pbmc.sh
```

### 5.3 Run BMMC RNA–ATAC

```bash
python preprocessing/scripts/preprocess_bmmc_rna_atac.py \
  --config preprocessing/configs/bmmc.yaml \
  2>&1 | tee preprocessing/logs/bmmc_rna_atac.log
```

### 5.4 Run BMMC RNA–Protein

```bash
python preprocessing/scripts/preprocess_bmmc_rna_protein.py \
  --config preprocessing/configs/bmmc.yaml \
  2>&1 | tee preprocessing/logs/bmmc_rna_protein.log
```

### 5.5 Run both BMMC pipelines sequentially

```bash
bash preprocessing/run_bmmc.sh
```

Logs are stored in:

```text
preprocessing/logs/
```

For example:

```bash
tail -f preprocessing/logs/pbmc_rna_atac.log
```

---

## 6. Generated data structure

The preprocessing scripts regenerate intermediate files from the raw public datasets rather than assuming that previously generated H5AD files already exist.

### 6.1 PBMC RNA–ATAC

Intermediate outputs:

```text
data/interim/PBMC/RNA_ATAC/
├── RNA_counts_qc.h5ad
├── ATAC_counts_qc.h5ad
├── ATAC_gas.h5ad
└── feature_aligned.h5ad
```

`ATAC_gas.h5ad` is the gene-level activity representation derived from peak-level ATAC data using the GTF annotation.

`feature_aligned.h5ad` retains genes that:

1. are highly variable in RNA or ATAC gene activity; and
2. are present in both RNA and ATAC gene-activity feature spaces.

Final split datasets:

```text
data/processed/PBMC/RNA_ATAC/results_ratio_loop/
├── summary_all_ratios.csv
├── single_000/
├── single_020/
├── single_040/
├── single_060/
├── single_080/
└── single_100/
```

Each successful ratio directory contains:

```text
train_rna_ref.h5ad
train_atac_full.h5ad
train_atac_activity.h5ad
val_query_atac.h5ad
val_true_rna.h5ad
val_atac_activity.h5ad
split_info.json
```

---

### 6.2 BMMC RNA–ATAC

Intermediate outputs:

```text
data/interim/BMMC/RNA_ATAC/
├── RNA_counts_qc.h5ad
├── ATAC_counts_qc.h5ad
├── ATAC_gas.h5ad
└── feature_aligned_rna_atac.h5ad
```

Final split datasets:

```text
data/processed/BMMC/RNA_ATAC/results_ratio_loop_rna_atac/
├── summary_all_ratios.csv
├── single_000/
├── single_020/
├── single_040/
├── single_060/
├── single_080/
└── single_100/
```

Each successful ratio directory contains:

```text
train_rna_ref.h5ad
train_atac_full.h5ad
train_atac_activity.h5ad
val_query_atac.h5ad
val_true_rna.h5ad
val_atac_activity.h5ad
split_info.json
```

---

### 6.3 BMMC RNA–Protein

Intermediate outputs:

```text
data/interim/BMMC/RNA_PROTEIN/
├── RNA_counts_qc.h5ad
├── protein_counts_qc.h5ad
└── feature_aligned_rna_protein.h5ad
```

Final split datasets:

```text
data/processed/BMMC/RNA_PROTEIN/results_ratio_loop_rna_protein/
├── summary_all_ratios.csv
├── single_000/
├── single_020/
├── single_040/
├── single_060/
├── single_080/
└── single_100/
```

Each successful ratio directory contains:

```text
train_rna_ref.h5ad
train_protein_ref.h5ad
train_protein_full.h5ad
val_query_protein.h5ad
val_true_rna.h5ad
val_true_protein.h5ad
split_info.json
```

---

## 7. Partial-pairing ratios

The preprocessing pipelines generate datasets with different fractions of single-modality cells:

```text
single_000    0% single-modality cells
single_020   20% single-modality cells
single_040   40% single-modality cells
single_060   60% single-modality cells
single_080   80% single-modality cells
single_100  100% single-modality cells
```

These datasets are used to evaluate ATLAS under different degrees of missing or unpaired multi-omics measurements.

The exact split assignments are stored in:

```text
split_info.json
```

and summary statistics are stored in:

```text
summary_all_ratios.csv
```

---

## 8. Reproducibility notes

- Raw public datasets should remain unchanged under `data/raw/`.
- External biological annotation resources should be placed under `data/external/annotations/`.
- Intermediate files generated by preprocessing belong under `data/interim/`.
- Final model-ready datasets belong under `data/processed/`.
- Preprocessing scripts should not depend on manually created intermediate H5AD files.
- Dataset-specific paths and preprocessing parameters are defined in YAML configuration files rather than being hard-coded in the Python scripts.
- Random seeds used for train/validation and partial-pairing splits are defined in the configuration files.
- The PBMC ATAC fragments file is kept for provenance/future analyses but is not consumed by the current PBMC preprocessing workflow.
- The BMMC HGNC/CD mapping file is optional and affects feature annotation only; it is not required for RNA–Protein normalization itself.

