# Input data

Place the following files in this directory (or point `SCRNA_DATA_DIR` at a
directory that contains them). Large RDS objects are not shipped with the code.

| File | Approx. size | Role |
|------|--------------|------|
| `seurat_object.rds` | ~600 MB | Preprocessed Seurat object (SCT-normalized, Harmony-corrected) |
| `clustassess_object.rds` | ~350 MB | ClustAssess stability object (20-cluster partition and UMAP) |
| `gene_background.csv` | ~95 KB | Custom g:Profiler background (included) |

If your deposited files use different names, either rename them to the defaults
above or set:

```bash
export SCRNA_SEURAT_RDS=your_seurat_file.rds
export SCRNA_CLUSTASSESS_RDS=your_clustassess_file.rds
```

Do not commit `*.rds` into version control (see the kit `.gitignore`).
