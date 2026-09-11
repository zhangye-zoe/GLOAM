# ATLAS Data Preprocessing

## 1. Download links and external files

### PBMC RNA–ATAC

Dataset: [10k Human PBMCs Multiome](https://www.10xgenomics.com/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-x-1-standard-2-0-0)

Download:

```text
pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
```

Save to:

```text
data/raw/PBMC/
```

### BMMC RNA–ATAC / RNA–Protein

Dataset: [GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122)

Download:

```text
GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz
GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz
```

Decompress:

```bash
gunzip GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz
gunzip GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz
```

Save to:

```text
data/raw/BMMC/
```

### External annotations

Required for PBMC/BMMC RNA–ATAC:

```text
data/external/annotations/gencode.v38.primary_assembly.annotation.gtf
```

Optional for BMMC RNA–Protein:

```text
data/external/annotations/hgnc_cd_gene_mapping.csv
```

The HGNC mapping file can be copied from the previous project:

```bash
cp /data5/zhangye/scMRDR/scripts/BMMC/CD_gene/group-471.csv \
   /data5/zhangye/ATLAS/data/external/annotations/hgnc_cd_gene_mapping.csv
```

---

## 2. Preprocessing scripts

```text
preprocessing/
├── configs/
│   ├── bmmc.yaml
│   └── pbmc.yaml
├── logs/
├── scripts/
│   ├── common.py
│   ├── preprocess_bmmc_rna_atac.py
│   ├── preprocess_bmmc_rna_protein.py
│   └── preprocess_pbmc_rna_atac.py
├── run_bmmc.sh
└── run_pbmc.sh
```

Run from the ATLAS root:

```bash
cd /data5/zhangye/ATLAS
```

PBMC RNA–ATAC:

```bash
bash preprocessing/run_pbmc.sh
```

BMMC RNA–ATAC and RNA–Protein:

```bash
bash preprocessing/run_bmmc.sh
```

Or run scripts individually:

```bash
python preprocessing/scripts/preprocess_pbmc_rna_atac.py \
  --config preprocessing/configs/pbmc.yaml

python preprocessing/scripts/preprocess_bmmc_rna_atac.py \
  --config preprocessing/configs/bmmc.yaml

python preprocessing/scripts/preprocess_bmmc_rna_protein.py \
  --config preprocessing/configs/bmmc.yaml
```

---

## 3. Generated directory structure

```text
data/
├── raw/
│   ├── PBMC/
│   │   └── pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
│   │   
│   └── BMMC/
│       ├── GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad
│       └── GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad
│
├── external/
│   └── annotations/
│       ├── gencode.v38.primary_assembly.annotation.gtf
│       └── hgnc_cd_gene_mapping.csv
│
├── interim/
│   ├── PBMC/
│   │   └── RNA_ATAC/
│   │       ├── RNA_counts_qc.h5ad
│   │       ├── ATAC_counts_qc.h5ad
│   │       ├── ATAC_gas.h5ad
│   │       └── feature_aligned.h5ad
│   │
│   └── BMMC/
│       ├── RNA_ATAC/
│       │   ├── RNA_counts_qc.h5ad
│       │   ├── ATAC_counts_qc.h5ad
│       │   ├── ATAC_gas.h5ad
│       │   └── feature_aligned_rna_atac.h5ad
│       │
│       └── RNA_PROTEIN/
│           ├── RNA_counts_qc.h5ad
│           ├── protein_counts_qc.h5ad
│           └── feature_aligned_rna_protein.h5ad
│
└── processed/
    ├── PBMC/
    │   └── RNA_ATAC/
    │       └── results_ratio_loop/
    │
    └── BMMC/
        ├── RNA_ATAC/
        │   └── results_ratio_loop_rna_atac/
        └── RNA_PROTEIN/
            └── results_ratio_loop_rna_protein/
```

Each `results_ratio_loop*` directory contains:

```text
single_000/
single_020/
single_040/
single_060/
single_080/
single_100/
summary_all_ratios.csv
```
