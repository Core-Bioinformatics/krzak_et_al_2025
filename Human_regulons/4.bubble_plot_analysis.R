library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

project_dir <- "/iss-scratch/CoreBioinformatics/rk720/human_sucnr1"

source(
  file.path(
    project_dir,
    "helper_functions.R"
  )
)

objects_dir <- file.path(
  project_dir,
  "objects",
  "seurat"
)

output_dir <- file.path(
  project_dir,
  "SUCNR1_bubble_plot_v3"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

genes_use <- c(
  "NFKB1",
  "NFKBIA",
  "HIF1A",
  "NLRP3",
  "CASP1",
  "ACOD1",
  "MRC1",
  "IL6",
  "IL10",
  "SPI1",
  "IRF8",
  "RUNX1",
  "P2RY12",
  "CX3CR1",
  "TREM2",
  "TYROBP",
  "ITGAM",
  "CSF1R",
  "C1QA",
  "TMEM119",
  "PARP14",
  "CTSS",
  "IFITM1",
  "MMP1",
  "BCL6",
  "CLEC7A",
  "CD68"
)

dataset_files <- c(
  Absinta = file.path(
    objects_dir,
    "absinta_microglia_obj.rds"
  ),
  Schirmer = file.path(
    objects_dir,
    "schirmer_microglia_obj.rds"
  ),
  MacNair = file.path(
    objects_dir,
    "macnair_microglia_obj.rds"
  )
)


calculate_dotplot_values <- function(
  obj,
  dataset_name,
  genes,
  positive_col = "SUCNR1_positive",
  homeostatic_markers = c(
    "P2RY12",
    "TMEM119"
  ),
  homeostatic_quantile = 0.75,
  n_resamples = 100,
  seed = 1234
) {
  if (!positive_col %in% colnames(obj@meta.data)) {
    stop(
      positive_col,
      " is missing from ",
      dataset_name
    )
  }

  if (!"RNA" %in% Assays(obj)) {
    stop("RNA assay missing from ", dataset_name)
  }

  obj <- NormalizeData(
    obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  counts <- get_assay_data_safe(
    obj = obj,
    assay = "RNA",
    layer = "counts"
  )

  normalized <- get_assay_data_safe(
    obj = obj,
    assay = "RNA",
    layer = "data"
  )

  missing_homeostatic_markers <- setdiff(
    homeostatic_markers,
    rownames(normalized)
  )

  if (length(missing_homeostatic_markers) > 0) {
    stop(
      dataset_name,
      ": missing homeostatic markers: ",
      paste(
        missing_homeostatic_markers,
        collapse = ", "
      )
    )
  }

  genes_present <- genes[
    genes %in% rownames(counts) &
      genes %in% rownames(normalized)
  ]

  genes_missing <- setdiff(
    genes,
    genes_present
  )

  if (length(genes_missing) > 0) {
    warning(
      dataset_name,
      ": missing genes: ",
      paste(
        genes_missing,
        collapse = ", "
      )
    )
  }

  cell_ids <- colnames(obj)

  sucnr1_positive <- as.logical(
    obj@meta.data[
      cell_ids,
      positive_col,
      drop = TRUE
    ]
  )

  positive_cells <- cell_ids[
    !is.na(sucnr1_positive) &
      sucnr1_positive
  ]

  negative_cells <- cell_ids[
    !is.na(sucnr1_positive) &
      !sucnr1_positive
  ]

  n_positive <- length(positive_cells)

  if (n_positive == 0) {
    stop(
      "No SUCNR1-positive cells in ",
      dataset_name
    )
  }

  if (length(negative_cells) < n_positive) {
    stop(
      dataset_name,
      ": fewer SUCNR1-negative cells than SUCNR1-positive cells."
    )
  }

  # Combined P2RY12/TMEM119 score among SUCNR1-negative cells.
  marker_expression <- t(
    as.matrix(
      normalized[
        homeostatic_markers,
        negative_cells,
        drop = FALSE
      ]
    )
  )

  # Scale each marker separately so neither marker dominates.
  marker_expression_scaled <- scale(
    marker_expression
  )

  marker_expression_scaled[
    !is.finite(marker_expression_scaled)
  ] <- 0

  homeostatic_score <- rowMeans(
    marker_expression_scaled
  )

  names(homeostatic_score) <- rownames(
    marker_expression
  )

  score_cutoff <- quantile(
    homeostatic_score,
    probs = homeostatic_quantile,
    na.rm = TRUE,
    names = FALSE
  )

  homeostatic_pool <- names(
    homeostatic_score
  )[
    homeostatic_score >= score_cutoff
  ]

  # Fallback if the top quartile is smaller than the
  # number of SUCNR1-positive cells.
  if (length(homeostatic_pool) < n_positive) {
    warning(
      dataset_name,
      ": top homeostatic quartile contains only ",
      length(homeostatic_pool),
      " cells. Using the top ",
      n_positive,
      " ranked cells instead."
    )

    homeostatic_pool <- names(
      sort(
        homeostatic_score,
        decreasing = TRUE
      )
    )[
      seq_len(n_positive)
    ]
  }

  message(
    dataset_name,
    ": ",
    n_positive,
    " SUCNR1+ cells; ",
    length(homeostatic_pool),
    " cells in the homeostatic control pool."
  )

  summarise_cells <- function(
    group_cells,
    group_name
  ) {
    percent_expression <- as.numeric(
      Matrix::rowMeans(
        counts[
          genes_present,
          group_cells,
          drop = FALSE
        ] > 0
      )
    ) * 100

    mean_expression <- as.numeric(
      Matrix::rowMeans(
        normalized[
          genes_present,
          group_cells,
          drop = FALSE
        ]
      )
    )

    data.frame(
      dataset = dataset_name,
      group = group_name,
      gene = genes_present,
      n_group_cells = length(group_cells),
      percent_expression = percent_expression,
      mean_expression = mean_expression,
      stringsAsFactors = FALSE
    )
  }

  # Fixed SUCNR1-positive group.
  positive_summary <- summarise_cells(
    group_cells = positive_cells,
    group_name = "SUCNR1+"
  )

  positive_summary$percent_expression_sd <- NA_real_
  positive_summary$mean_expression_sd <- NA_real_
  positive_summary$homeostatic_pool_size <-
    length(homeostatic_pool)
  positive_summary$n_resamples <- 1L

  # Repeated equal-size sampling from the homeostatic pool.
  set.seed(seed)

  homeostatic_resamples <- bind_rows(
    lapply(
      seq_len(n_resamples),
      function(resample_id) {
        sampled_cells <- sample(
          homeostatic_pool,
          size = n_positive,
          replace = FALSE
        )

        summary <- summarise_cells(
          group_cells = sampled_cells,
          group_name = "Homeostatic"
        )

        summary$resample_id <- resample_id

        summary
      }
    )
  )

  homeostatic_summary <- homeostatic_resamples %>%
  group_by(
    dataset,
    group,
    gene
  ) %>%
  summarise(
    percent_expression_mean = mean(
      percent_expression,
      na.rm = TRUE
    ),
    percent_expression_sd = sd(
      percent_expression,
      na.rm = TRUE
    ),
    mean_expression_mean = mean(
      mean_expression,
      na.rm = TRUE
    ),
    mean_expression_sd = sd(
      mean_expression,
      na.rm = TRUE
    ),
    n_group_cells = first(n_group_cells),
    .groups = "drop"
  ) %>%
  mutate(
    percent_expression = percent_expression_mean,
    mean_expression = mean_expression_mean,
    homeostatic_pool_size = length(homeostatic_pool),
    n_resamples = n_resamples
  ) %>%
  select(
    -percent_expression_mean,
    -mean_expression_mean
  )

  bind_rows(
    homeostatic_summary,
    positive_summary
  )
}

dotplot_data <- bind_rows(
  lapply(
    names(dataset_files),
    function(dataset_name) {
      message("Processing ", dataset_name)

      obj <- readRDS(
        dataset_files[[dataset_name]]
      )

      calculate_dotplot_values(
        obj = obj,
        dataset_name = dataset_name,
        genes = genes_use
      )
    }
  )
)


# Scale the mean expression separately for each gene across the six dataset/group combinations.
dotplot_data <- dotplot_data %>%
  group_by(gene) %>%
  mutate(
    relative_expression = {
      gene_sd <- sd(
        mean_expression,
        na.rm = TRUE
      )

      if (is.na(gene_sd) || gene_sd == 0) {
        rep(0, n())
      } else {
        (
          mean_expression -
            mean(
              mean_expression,
              na.rm = TRUE
            )
        ) / gene_sd
      }
    }
  ) %>%
  ungroup()

# Clip extreme values only for colour visualization.
dotplot_data$relative_expression_plot <- pmax(
  -2,
  pmin(
    2,
    dotplot_data$relative_expression
  )
)

dotplot_data$dataset <- factor(
  dotplot_data$dataset,
  levels = c(
    "Absinta",
    "Schirmer",
    "MacNair"
  )
)

dotplot_data$group <- factor(
  dotplot_data$group,
  levels = c(
    "Homeostatic",
    "SUCNR1+"
  )
)

dotplot_data <- dotplot_data %>%
  mutate(
    group_label = paste0(
      as.character(group),
      "\n(n=",
      n_group_cells,
      ")"
    )
  )

present_gene_order <- genes_use[
  genes_use %in% unique(dotplot_data$gene)
]

# Reverse so the first gene appears at the top.
dotplot_data$gene <- factor(
  dotplot_data$gene,
  levels = rev(present_gene_order)
)


bubble_plot <- ggplot(
  dotplot_data,
  aes(
    x = group_label,
    y = gene
  )
) +
  geom_point(
    aes(
      size = percent_expression,
      color = relative_expression_plot
    )
  ) +
  facet_grid(
    cols = vars(dataset),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_color_gradient2(
    low = "royalblue4",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-2, 2),
    name = "Relative\nexpression"
  ) +
  scale_size_area(
    max_size = 8,
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = "Percent\nexpression"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "italic",
      size = 12
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.spacing.x = unit(
      1,
      "lines"
    ),
    plot.subtitle = element_text(
      size = 10,
      lineheight = 1.1,
      margin = margin(b = 10)
    )
  ) +
  labs(
  title = paste0(
    "Gene expression in matched homeostatic and ",
    "SUCNR1+ microglia"
  ),
  subtitle = paste0(
    "Homeostatic controls: equal-size samples from the top 25% of SUCNR1− cells\n",
    "ranked by combined P2RY12/TMEM119 expression (100 resamples).\n",
    "Dot size: % detected; colour: mean log-normalized expression scaled separately per gene.\n",
    "Compare colours horizontally within each gene row."
  )
)


ggsave(
  filename = file.path(
    output_dir,
    "SUCNR1_positive_negative_gene_bubble_plot.png"
  ),
  plot = bubble_plot,
  width = 11,
  height = 12,
  dpi = 300
)

ggsave(
  filename = file.path(
    output_dir,
    "SUCNR1_positive_negative_gene_bubble_plot.pdf"
  ),
  plot = bubble_plot,
  width = 11,
  height = 12
)

write.csv(
  dotplot_data,
  file.path(
    output_dir,
    "SUCNR1_positive_negative_gene_bubble_plot_values.csv"
  ),
  row.names = FALSE
)