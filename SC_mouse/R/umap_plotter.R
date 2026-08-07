# =============================================================================
# UMAP Plotting Utilities for Seurat Objects
# =============================================================================
#
# Usage:
#   source("umap_plotter.R")
#   plot_umap_metadata(seurat_obj, feature = "cluster_col")
#
# Requirements:
#   - ggplot2
#   - colorspace (for make_umap_palette)
#
# =============================================================================

extract_umap_limits <- function(seurat_obj,
                                embedding_name = "umap",
                                padding_frac = 0.03) {
  stopifnot(padding_frac >= 0)

  umap_mat <- as.data.frame(seurat_obj@reductions[[embedding_name]]@cell.embeddings)
  main_x <- colnames(umap_mat)[1]
  main_y <- colnames(umap_mat)[2]

  add_padding <- function(x) {
    xr <- range(x, na.rm = TRUE)
    span <- diff(xr)
    if (!is.finite(span) || span == 0) span <- 1
    pad <- span * padding_frac
    c(xr[1] - pad, xr[2] + pad)
  }

  list(
    x = add_padding(umap_mat[[main_x]]),
    y = add_padding(umap_mat[[main_y]])
  )
}

save_plot_with_fixed_legend <- function(plot_obj,
                                        file,
                                        width,
                                        height,
                                        legend_width_cm = 5,
                                        legend_gap_cm = 0.4) {
  plot_no_legend <- plot_obj + ggplot2::theme(legend.position = "none")
  plot_grob <- ggplot2::ggplotGrob(plot_no_legend)
  full_grob <- ggplot2::ggplotGrob(plot_obj)

  guide_idx <- which(full_grob$layout$name == "guide-box-right")
  if (length(guide_idx) == 0) {
    guide_idx <- which(grepl("guide-box", full_grob$layout$name, fixed = TRUE))
  }

  grDevices::pdf(
    file = file,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()

  if (length(guide_idx) == 0) {
    grid::grid.draw(plot_grob)
    return(invisible(file))
  }

  legend_grob <- full_grob$grobs[[guide_idx[1]]]
  layout <- grid::grid.layout(
    nrow = 1,
    ncol = 3,
    widths = grid::unit.c(
      grid::unit(1, "null"),
      grid::unit(legend_gap_cm, "cm"),
      grid::unit(legend_width_cm, "cm")
    )
  )

  grid::pushViewport(grid::viewport(layout = layout))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.draw(plot_grob)
  grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
  grid::grid.draw(legend_grob)
  grid::popViewport(2)

  invisible(file)
}

sort_feature_levels <- function(feature_levels) {
  nums <- suppressWarnings(as.numeric(feature_levels))
  if (!any(is.na(nums))) {
    feature_levels[order(nums)]
  } else {
    sort(feature_levels)
  }
}

# Default palette for clusters (up to 21 types, hand-picked for aesthetics)
umap21_fixed <- c(
  "#CBD5CD", # sage grey
  "#EEA2AD", # dusty pink
  "#607B8B", # blue-grey
  "#8B7B8B", # muted purple-grey
  "#4682B4", # steel blue
  "#CDC8B1", # warm grey-beige    # split one
  "#BCD2EE", # powder blue
  "#79CDCD", # turquoise
  "#5F9EA0", # cadet teal
  "#698B69", # muted green
  "#C1E1C1", # pale mint
  "#7FA58A", # soft green-teal
  "#CDAA7D", # tan
  "#F4C27D", # warm sand
  "#D9A441", # muted mustard
  "#CDB5CD", # lavender
  "#123524", # phthalo green
  "#A67C52", # muted caramel
  "#8B4513", # deep brown
  "#B48EAD", # mauve
  "#D0A9A9"  # dusty rose
)

#' Generate a perceptually-optimized qualitative palette for UMAP plots
#'
#' Uses HCL color space for perceptual uniformity, then reorders colors via
#' farthest-point sampling in CIE Lab space to maximize sequential distinguishability.
#'
#' @param k Number of colors to generate (must be >= 2)
#' @param l Lightness (0-100). Higher = more pastel. Default 72.
#' @param c Chroma/saturation (0-100+). Higher = more vivid. Default 45.
#' @param h_start Starting hue (0-360). Default 15 (warm orange-ish start).
#' @param seed Random seed for reproducible starting point in farthest-point sampling.
#' @return Character vector of k hex colors, ordered for maximum sequential contrast.
#'
#' @details Requires the colorspace package. The algorithm:
#'   1. Generate k evenly-spaced colors in HCL space
#'   2. Convert to CIE Lab (perceptually uniform)
#'   3. Reorder via farthest-point sampling: each successive color is chosen
#'      to be maximally distant from all previously selected colors
#'
#' @examples
#' make_umap_palette(10)
#' make_umap_palette(20, l = 65, c = 55)  # darker, more saturated
#'
make_umap_palette <- function(k,
                              l = 72,
                              c = 45,
                              h_start = 15,
                              seed = 1) {
  stopifnot(k >= 2)
  stopifnot(l >= 0 && l <= 100)
  stopifnot(c >= 0)
  stopifnot(h_start >= 0 && h_start <= 360)

  if (!requireNamespace("colorspace", quietly = TRUE)) {
    stop("Package 'colorspace' is required for make_umap_palette(). Install with: install.packages('colorspace')")
  }

  set.seed(seed)

  # Generate k evenly-spaced colors in HCL space
  base <- colorspace::qualitative_hcl(
    n = k,
    h = c(h_start, h_start + 360),
    c = c,
    l = l
  )

  # Convert to Lab for perceptually-uniform distance calculations
  rgb_mat <- grDevices::col2rgb(base) / 255
  lab_mat <- grDevices::convertColor(t(rgb_mat), from = "sRGB", to = "Lab")

  # Farthest-point sampling: reorder colors to maximize sequential separation
  # This ensures adjacent colors in the legend are maximally distinct
  idx <- integer(k)
  idx[1] <- sample.int(k, 1)
  remaining <- setdiff(seq_len(k), idx[1])

  for (i in 2:k) {
    # For each remaining color, compute min distance to any selected color
    selected_lab <- lab_mat[idx[seq_len(i - 1)], , drop = FALSE]

    min_dists <- vapply(remaining, function(j) {
      # Euclidean distance in Lab space to each selected color
      diffs <- sweep(selected_lab, 2, lab_mat[j, ])
      min(sqrt(rowSums(diffs^2)))
    }, numeric(1))

    # Select the color with maximum min-distance (farthest from all selected)
    best_idx <- which.max(min_dists)
    idx[i] <- remaining[best_idx]
    remaining <- remaining[-best_idx]
  }

  base[idx]
}

#' Plot UMAP with customizable coloring and filtering
#'
#' @param seurat_obj    Seurat object
#' @param embedding_name Name of reduction (default 'umap')
#' @param feature Metadata column to color by
#' @param filters List of filters for metadata columns
#' @param palette_colors Optional vector of colors (overrides palette_mode)
#' @param palette_mode How to generate colors when palette_colors is NULL:
#'   - "fixed": use hand-picked umap21_fixed palette (up to 21 colors)
#'   - "auto": dynamically generate colors via make_umap_palette() (any k)
#' @param palette_l Lightness for auto palette (0-100, higher = more pastel)
#' @param palette_c Chroma for auto palette (higher = more saturated)
#' @param add_title Add the plot title?
#' @param plot_title Custom plot title (NULL = auto-generate)
#' @param legend_title Legend title for the color scale (NULL = use `feature`)
#' @param append_filters_to_title If filters are applied, append them to the title?
#' @param pdf_file_out Optional PDF output path. If provided, the plot is saved as a PDF.
#' @param pdf_width PDF width in inches (default 20)
#' @param pdf_height PDF height in inches (default 20)
#' @param pdf_autoscale_style_sizes If TRUE and pdf_file_out is set, autoscale sizes for
#'   larger-format export (unless you explicitly pass size arguments yourself).
#' @param pdf_autoscale_base_inches Baseline inches for the provided style defaults (default 15)
#' @param point_size Size of scatter plot points
#' @param alpha Alpha of foreground points
#' @param background_alpha Background points alpha
#' @param background_color Background color
#' @param base_text_size Base text size passed to ggplot2::theme_classic()
#' @param legend_title_size Legend title text size (NULL = ggplot default)
#' @param legend_text_size Legend item text size (NULL = ggplot default)
#' @param axis_title_size Axis title text size (NULL = ggplot default)
#' @param axis_text_size Axis tick label text size (NULL = ggplot default)
#' @param plot_title_size Plot title text size (NULL = ggplot default)
#' @param legend_key_point_size Increase legend point size (NULL = default)
#' @param fixed_limits Optional list with numeric vectors `x` and `y` for shared axis limits
#' @param coord_expand Passed to ggplot2::coord_fixed(expand = ...)
#' @param legend_width_cm Reserved width for the legend column in exported PDFs
#' @param legend_gap_cm Gap between panel and legend in exported PDFs
#' @param order_method How to order points for plotting occlusion:
#'   - "density": per-point density ordering (dense first, sparse last)
#'   - "shuffle": random shuffle (legacy)
#'   - "none": original row order
#' @param density_bins Grid resolution for density ordering
#' @param seed Seed for shuffle/palette generation
#'
plot_umap_metadata <- function(seurat_obj,
                              embedding_name = "umap",
                              feature = "stable_20_clusters",
                              filters = list(),
                              palette_colors = NULL,
                              palette_mode = c("fixed", "auto"),
                              palette_l = 72,
                              palette_c = 45,
                              add_title = TRUE,
                              plot_title = NULL,
                              legend_title = NULL,
                              append_filters_to_title = NULL,
                              pdf_file_out = NULL,
                              pdf_width = 20,
                              pdf_height = 20,
                              pdf_autoscale_style_sizes = TRUE,
                              pdf_autoscale_base_inches = 15,
                              point_size = 0.3,
                              alpha = 0.6,
                              background_alpha = 0.2,
                              background_color = "#BDBDBD",
                              base_text_size = 11,
                              legend_title_size = 14,
                              legend_text_size = 12,
                              axis_title_size = 14,
                              axis_text_size = 12,
                              plot_title_size = 14,
                              legend_key_point_size = 4,
                              fixed_limits = NULL,
                              coord_expand = FALSE,
                              legend_width_cm = 5,
                              legend_gap_cm = 0.4,
                              order_method = c("density", "shuffle", "none"),
                              density_bins = 100,
                              seed = 42) {

  # This line checks that 'order_method' matches one of the allowed options
  # ("density", "shuffle", or "none") and sets 'order_method' to the matched value.
  # If 'order_method' is not specified, "density" will be used as the default.
  order_method <- match.arg(order_method)
  palette_mode <- match.arg(palette_mode)

  if (!is.null(pdf_file_out) && pdf_autoscale_style_sizes) {
    scale_factor <- min(pdf_width, pdf_height) / pdf_autoscale_base_inches

    # Baseline values were tuned for a 15x15-inch PDF.
    if (missing(point_size)) {
      point_size <- 0.8 * scale_factor
    }
    if (missing(legend_text_size)) {
      legend_text_size <- 16 * scale_factor
    }
    if (missing(legend_key_point_size)) {
      legend_key_point_size <- 7 * scale_factor
    }
    if (missing(legend_title_size)) {
      legend_title_size <- 18 * scale_factor
    }
    if (missing(plot_title_size)) {
      plot_title_size <- 18 * scale_factor
    }

    # Give the overall theme a sensible bump unless the caller set it.
    if (missing(base_text_size)) {
      base_text_size <- 14 * scale_factor
    }
  }

  # Extract UMAP embedding
  umap_mat <- as.data.frame(seurat_obj@reductions[[embedding_name]]@cell.embeddings)
  # Get metadata columns needed for plotting and filtering
  meta_df <- seurat_obj[[]]

  # Merge embedding and metadata
  umap_df <- cbind(umap_mat, meta_df)

  # Defensive: check feature exists
  if (!(feature %in% colnames(umap_df))) {
    stop(paste("Feature", feature, "not found in cell metadata."))
  }

  main_x <- colnames(umap_mat)[1]
  main_y <- colnames(umap_mat)[2]

  if (!is.null(fixed_limits)) {
    stopifnot(is.list(fixed_limits))
    stopifnot(all(c("x", "y") %in% names(fixed_limits)))
    stopifnot(length(fixed_limits$x) == 2, length(fixed_limits$y) == 2)
  }

  filters_applied <- NULL

  # Identify which cells pass the filters, if any
  if (length(filters) > 0) {
    filt_idx <- rep(TRUE, nrow(umap_df))
    filters_strings <- character(length(filters))
    idx <- 1
    for (filt_name in names(filters)) {
      filt_vals <- filters[[filt_name]]
      filt_idx <- filt_idx & (umap_df[[filt_name]] %in% filt_vals)
      filters_strings[idx] <- paste0(filt_name, "=", paste(filt_vals, collapse = "|"))
      idx <- idx + 1
    }
    filters_applied <- paste(filters_strings, collapse = ", ")
  } else {
    filt_idx <- rep(TRUE, nrow(umap_df))
  }

  umap_df$highlighted <- filt_idx # logical: TRUE=pass filter, FALSE=background

  # Split into background and highlighted data frames
  umap_df_bg <- umap_df[!umap_df$highlighted, , drop = FALSE]
  umap_df_hl <- umap_df[umap_df$highlighted, , drop = FALSE]

  # -----------------------------------------------------------------
  # Point ordering for occlusion management
  # -----------------------------------------------------------------
  if (nrow(umap_df_hl) > 0) {
    if (order_method == "density") {
      # Per-point density ordering: dense regions first, sparse last
      # This exposes rare/sparse points by drawing them on top

      xr <- range(umap_df_hl[[main_x]], na.rm = TRUE)
      yr <- range(umap_df_hl[[main_y]], na.rm = TRUE)

      # Bin each point into a grid cell
      bx <- pmin(density_bins, pmax(1, floor((umap_df_hl[[main_x]] - xr[1]) / diff(xr) * density_bins) + 1))
      by <- pmin(density_bins, pmax(1, floor((umap_df_hl[[main_y]] - yr[1]) / diff(yr) * density_bins) + 1))
      bin_id <- paste(bx, by, sep = "_")

      # Count points per bin (local density)
      bin_counts <- table(bin_id)
      umap_df_hl$.density <- as.integer(bin_counts[bin_id])

      # Compute cluster/feature sizes for tie-breaking (small clusters on top)
      feature_vals_hl <- umap_df_hl[[feature]]
      if (!is.numeric(feature_vals_hl)) {
        feature_sizes <- table(feature_vals_hl)
        umap_df_hl$.feature_size <- as.integer(feature_sizes[as.character(feature_vals_hl)])
      } else {
        # For numeric features, no cluster size tie-breaker
        umap_df_hl$.feature_size <- 0
      }

      # Sort: dense first (-density), big clusters first (-feature_size)
      # Result: sparse points from small clusters end up on top (drawn last)
      umap_df_hl <- umap_df_hl[order(-umap_df_hl$.density, -umap_df_hl$.feature_size), , drop = FALSE]

      # Clean up temporary columns
      umap_df_hl$.density <- NULL
      umap_df_hl$.feature_size <- NULL

    } else if (order_method == "shuffle") {
      # Legacy: random shuffle to remove systematic bias
      set.seed(seed)
      umap_df_hl <- umap_df_hl[sample(nrow(umap_df_hl)), , drop = FALSE]
    }
    # order_method == "none": keep original row order
  }

  # Prep palette
  feature_vals <- umap_df[[feature]]
  feature_vals_highlight <- umap_df_hl[[feature]]

  if (!is.null(palette_colors)) {
    # User-provided palette takes precedence
    color_scale <- ggplot2::scale_color_manual(values = palette_colors)
  } else if (is.numeric(feature_vals_highlight)) {
    # Continuous feature: viridis
    color_scale <- ggplot2::scale_color_viridis_c()
  } else {
    # Categorical feature: choose based on palette_mode
    feature_levels <- unique(as.character(feature_vals_highlight))
    n_levels <- length(feature_levels)
    sorted_levels <- sort_feature_levels(feature_levels)

    if (palette_mode == "fixed") {
      # Use hand-picked umap21_fixed palette
      if (n_levels <= length(umap21_fixed)) {
        manual_colors <- umap21_fixed[seq_len(n_levels)]
      } else {
        # Fallback for > 21 levels: blue-yellow-red interpolation
        blue_yellow_red <- c("#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8",
                             "#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026")
        manual_colors <- colorRampPalette(blue_yellow_red)(n_levels)
      }
    } else {
      # palette_mode == "auto": dynamically generate via make_umap_palette
      if (n_levels >= 2) {
        manual_colors <- make_umap_palette(n_levels, l = palette_l, c = palette_c, seed = seed)
      } else {
        # Edge case: single level
        manual_colors <- "#4682B4"
      }
    }

    names(manual_colors) <- sorted_levels
    color_scale <- ggplot2::scale_color_manual(values = manual_colors)
  }

  # Build title string
  if (is.null(append_filters_to_title)) {
    append_filters_to_title <- is.null(plot_title)
  }

  title_str <- if (is.null(plot_title)) paste0("UMAP: ", feature) else plot_title
  if (!is.null(filters_applied) && add_title && append_filters_to_title) {
    title_str <- paste0(title_str, " | Filters: ", filters_applied)
  }

  legend_title_str <- if (is.null(legend_title)) feature else legend_title

  # Prep plot (background first, then ordered highlight points)
  p <- ggplot2::ggplot() +
    # draw background (all cells not passing the filter, or just all if no filter)
    ggplot2::geom_point(
      data = umap_df_bg,
      ggplot2::aes(x = .data[[main_x]], y = .data[[main_y]]),
      color = background_color,
      size = point_size,
      alpha = background_alpha,
      na.rm = TRUE
    ) +
    # draw filtered points on top, colored by feature (density-ordered)
    ggplot2::geom_point(
      data = umap_df_hl,
      ggplot2::aes(x = .data[[main_x]], y = .data[[main_y]], color = .data[[feature]]),
      size = point_size,
      alpha = alpha,
      na.rm = TRUE
    ) +
    color_scale +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = if (is.null(fixed_limits)) NULL else fixed_limits$x,
      ylim = if (is.null(fixed_limits)) NULL else fixed_limits$y,
      expand = coord_expand
    ) +
    ggplot2::theme_classic(base_size = base_text_size) +
    ggplot2::theme(
      plot.title = if (is.null(plot_title_size)) {
        ggplot2::element_text(hjust = 0.5)
      } else {
        ggplot2::element_text(hjust = 0.5, size = plot_title_size)
      },
      legend.position = "right"
    )

  if (!is.null(legend_title_size)) {
    p <- p + ggplot2::theme(legend.title = ggplot2::element_text(size = legend_title_size))
  }

  if (!is.null(legend_text_size)) {
    p <- p + ggplot2::theme(legend.text = ggplot2::element_text(size = legend_text_size))
  }

  if (!is.null(axis_title_size)) {
    p <- p + ggplot2::theme(axis.title = ggplot2::element_text(size = axis_title_size))
  }

  if (!is.null(axis_text_size)) {
    p <- p + ggplot2::theme(axis.text = ggplot2::element_text(size = axis_text_size))
  }

  if (!is.null(legend_key_point_size) && !is.numeric(feature_vals_highlight)) {
    p <- p + ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = legend_key_point_size, alpha = 1))
    )
  }

  if (add_title) {
    p <- p + ggplot2::labs(title = title_str, color = legend_title_str)
  }

  if (!is.null(pdf_file_out)) {
    save_plot_with_fixed_legend(
      plot_obj = p,
      file = pdf_file_out,
      width = pdf_width,
      height = pdf_height,
      legend_width_cm = legend_width_cm,
      legend_gap_cm = legend_gap_cm
    )
  }

  return(p)
}
