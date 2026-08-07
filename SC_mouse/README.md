# scRNA-seq analysis scripts (Krzak et al.)

Scripts used for the Online Methods analyses of the ex vivo mouse spinal-cord
single-cell RNA-seq dataset (Sucnr1; Ctrl and C-EAE).

This kit starts from preprocessed Seurat and ClustAssess objects (see `data/`).
It does not rebuild those objects from raw counts.

## Layout

```
SC_mouse/
  environment.yml   # conda environment
  install_extra.R   # ClustAssess (GitHub)
  README.md
  R/                # analysis scripts
  data/             # inputs (large RDS not bundled)
  output/           # default results directory
```

## Install

```bash
cd SC_mouse
conda env create -f environment.yml
conda activate krzak_scrna_methods
Rscript install_extra.R
```

Package versions used when preparing these scripts:

| Software | Version |
|----------|---------|
| R | 4.4.3 |
| Seurat | 5.4.0 |
| ClustAssess | 1.1.1 (GitHub ref `release-1.2.0`) |
| harmony | 1.2.4 |
| gprofiler2 | 0.2.4 |
| dittoSeq | 1.18.0 |

## Data

See [`data/README.md`](data/README.md). Place in `data/`:

- `seurat_object.rds`
- `clustassess_object.rds`
- `gene_background.csv` (included; g:Profiler background)

Optional overrides:

```bash
export SCRNA_METHODS_ROOT=/path/to/SC_mouse
export SCRNA_DATA_DIR=/path/to/inputs
export SCRNA_OUT_DIR=/path/to/outputs
export SCRNA_SEURAT_RDS=seurat_object.rds
export SCRNA_CLUSTASSESS_RDS=clustassess_object.rds
```

## Run order

From the `SC_mouse/` directory:

```bash
Rscript R/export_umaps.R

Rscript R/cell_abundance_fisher.R
Rscript R/cell_abundance_stats.R

Rscript R/deg_cluster_markers.R
Rscript R/heatmap_cluster_markers.R

Rscript R/deg_ctrl_vs_ko.R
Rscript R/go_ctrl_vs_ko.R          # needs network access (g:Profiler API)

Rscript R/eae_bubble_plots_global_microglia.R
Rscript R/eae_bubble_plots_clusters.R

# Optional
Rscript R/arg1_cluster_expression.R
```

Shared helpers (sourced by the scripts above):

- `R/00_paths.R` — paths and output helpers
- `R/load_annotated_object.R` — load Seurat + ClustAssess partition + marker voting
- `R/marker_genes.R` — `add_celltype_markers()`
- `R/umap_plotter.R` — UMAP plotting helpers

## Notes

- Differential-expression scripts are the longest-running step and need substantial RAM (prefer ≥64 GB).
- Results are written under `output/<analysis_name>/`.
