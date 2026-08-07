# Install GitHub-only R packages not provided by environment.yml.
# Run after: conda env create -f environment.yml && conda activate krzak_scrna_methods
#
# ClustAssess is required by the shared loader and UMAP export.

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

ref <- "release-1.2.0"
message("Installing ClustAssess from GitHub @", ref, " ...")
remotes::install_github(
  paste0("Core-Bioinformatics/ClustAssess@", ref),
  upgrade = "never"
)

message("ClustAssess: ", as.character(utils::packageVersion("ClustAssess")))
