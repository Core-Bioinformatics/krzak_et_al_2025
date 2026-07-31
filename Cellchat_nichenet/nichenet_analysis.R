library(AnnotationDbi)
library(GO.db)
library(org.Mm.eg.db)
library(nichenetr)
library(Seurat)
library(tidyverse)

object_file <- "/iss-scratch/CoreBioinformatics/rk720/nichenet_stuff/sucnr1_nichenet_object.rds"
network_dir <- "/iss-scratch/CoreBioinformatics/rk720/nichenet_stuff/useful_rds"
output_dir <- "/iss-scratch/CoreBioinformatics/rk720/nichenet_stuff/mechanistic_nichenet_corrected"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

conditions <- c("EAE_ctrl", "EAE_ko")
sender_celltypes_requested <- c("Microglia_Homeostatic", "Microglia_DAM")

minimum_cells_per_condition_primary <- 20
minimum_cells_per_condition_exploratory <- 2
exploratory_receivers_requested <- c("Neutrophils")

min_pct_expr <- 0.05
interest_min_pct_expr <- 0.01

interest_ligands <- c("Cxcl2")
interest_receptors <- c("Cxcr2", "Ackr1")

lfc_cutoff <- 0.25
padj_cutoff <- 0.05
minimum_geneset_size <- 5

# load data and normalize
seurat_obj <- readRDS(object_file)

stopifnot(all(c("celltype_markers", "eae_condition") %in% colnames(seurat_obj@meta.data)))
stopifnot("RNA" %in% Assays(seurat_obj))

DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "celltype_markers"

keep_cells <- rownames(seurat_obj@meta.data)[
  seurat_obj$eae_condition %in% conditions &
  !is.na(seurat_obj$celltype_markers) &
  as.character(seurat_obj$celltype_markers) != "Unassigned"
]

eae_obj <- subset(seurat_obj, cells = keep_cells)

if (inherits(eae_obj[["RNA"]], "Assay5")) {
  eae_obj <- JoinLayers(eae_obj, assay = "RNA")
}

eae_obj <- NormalizeData(eae_obj, assay = "RNA", verbose = FALSE)

DefaultAssay(eae_obj) <- "RNA"

eae_obj$celltype_markers <- droplevels(factor(eae_obj$celltype_markers))
eae_obj$eae_condition <- factor(eae_obj$eae_condition, levels = conditions)

Idents(eae_obj) <- "celltype_markers"


# Determine valid Senders and Receivers
cell_counts <- eae_obj@meta.data %>%
  dplyr::count(celltype_markers, eae_condition, name = "n_cells") %>%
  tidyr::complete(celltype_markers, eae_condition = conditions, fill = list(n_cells = 0)) %>%
  tidyr::pivot_wider(names_from = eae_condition, values_from = n_cells, values_fill = 0) %>%
  dplyr::arrange(celltype_markers)

print(cell_counts)

primary_valid_celltypes <- cell_counts %>%
  dplyr::filter(
    EAE_ctrl >= minimum_cells_per_condition_primary,
    EAE_ko >= minimum_cells_per_condition_primary
  ) %>%
  dplyr::pull(celltype_markers) %>%
  as.character()

exploratory_valid_celltypes <- cell_counts %>%
  dplyr::filter(
    EAE_ctrl >= minimum_cells_per_condition_exploratory,
    EAE_ko >= minimum_cells_per_condition_exploratory
  ) %>%
  dplyr::pull(celltype_markers) %>%
  as.character()

sender_celltypes <- intersect(sender_celltypes_requested, primary_valid_celltypes)

primary_receiver_celltypes <- setdiff(
  primary_valid_celltypes,
  sender_celltypes_requested
)

exploratory_receiver_celltypes <- intersect(
  exploratory_receivers_requested,
  exploratory_valid_celltypes
)

exploratory_receiver_celltypes <- setdiff(
  exploratory_receiver_celltypes,
  primary_receiver_celltypes
)

receiver_celltypes <- union(
  primary_receiver_celltypes,
  exploratory_receiver_celltypes
)

if (length(sender_celltypes) == 0) {
  stop("No requested senders meet the primary minimum cell-count threshold.")
}

if (length(receiver_celltypes) == 0) {
  stop("No receiver populations meet the minimum cell-count thresholds.")
}

receiver_classification <- tibble(receiver = receiver_celltypes) %>%
  dplyr::mutate(
    analysis_class = ifelse(
      receiver %in% primary_receiver_celltypes,
      "Primary",
      "Exploratory"
    )
  )

message("Valid Senders: ", paste(sender_celltypes, collapse = ", "))
message("Primary Receivers: ", paste(primary_receiver_celltypes, collapse = ", "))
message("Exploratory Receivers: ", paste(exploratory_receiver_celltypes, collapse = ", "))

# Load NicheNet Priors
ligand_target_matrix <- readRDS(
  file.path(network_dir, "ligand_target_matrix.rds")
)

lr_network <- readRDS(
  file.path(network_dir, "lr_network.rds")
) %>%
  dplyr::distinct(from, to)

nichenet_ligand_universe <- intersect(
  unique(lr_network$from),
  colnames(ligand_target_matrix)
)

nichenet_target_universe <- rownames(ligand_target_matrix)


# Define GO Cytokine and Chemokine Ligands
cytokine_go_root <- "GO:0005125"

cytokine_go_terms <- unique(
  c(
    cytokine_go_root,
    as.list(GO.db::GOMFOFFSPRING)[[cytokine_go_root]]
  )
)

cytokine_go_terms <- cytokine_go_terms[!is.na(cytokine_go_terms)]

cytokine_annotations <- suppressMessages(
  AnnotationDbi::select(
    x = org.Mm.eg.db,
    keys = cytokine_go_terms,
    keytype = "GOALL",
    columns = c("SYMBOL", "ONTOLOGYALL")
  )
) %>%
  tibble::as_tibble() %>%
  dplyr::filter(ONTOLOGYALL == "MF", !is.na(SYMBOL))

cytokine_ligands <- cytokine_annotations %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  intersect(nichenet_ligand_universe) %>%
  sort()

message("GO-defined cytokine/chemokine ligands: ", length(cytokine_ligands))


# Helper Functions
get_cells <- function(object, celltype, condition = NULL) {
  keep <- object$celltype_markers == celltype

  if (!is.null(condition)) {
    keep <- keep & object$eae_condition == condition
  }

  colnames(object)[keep]
}

get_expression_fraction <- function(object, celltype, condition = NULL) {
  cells <- get_cells(
    object = object,
    celltype = celltype,
    condition = condition
  )

  if (length(cells) == 0) {
    return(numeric())
  }

  expression_matrix <- GetAssayData(
    object,
    assay = "RNA",
    layer = "data"
  )[, cells, drop = FALSE]

  Matrix::rowMeans(expression_matrix > 0)
}

get_expressed_genes <- function(
  object,
  celltype,
  condition = NULL,
  expression_threshold = 0.05
) {
  expression_fraction <- get_expression_fraction(
    object = object,
    celltype = celltype,
    condition = condition
  )

  names(
    expression_fraction[
      expression_fraction >= expression_threshold
    ]
  )
}

get_eligible_sender_genes <- function(object, celltype, condition) {
  expression_fraction <- get_expression_fraction(
    object = object,
    celltype = celltype,
    condition = condition
  )

  standard_genes <- names(
    expression_fraction[
      expression_fraction >= min_pct_expr
    ]
  )

  interest_genes_passing <- intersect(
    interest_ligands,
    names(
      expression_fraction[
        expression_fraction >= interest_min_pct_expr
      ]
    )
  )

  union(standard_genes, interest_genes_passing)
}

get_eligible_receiver_genes <- function(object, celltype, condition) {
  expression_fraction <- get_expression_fraction(
    object = object,
    celltype = celltype,
    condition = condition
  )

  standard_genes <- names(
    expression_fraction[
      expression_fraction >= min_pct_expr
    ]
  )

  interest_genes_passing <- intersect(
    interest_receptors,
    names(
      expression_fraction[
        expression_fraction >= interest_min_pct_expr
      ]
    )
  )

  union(standard_genes, interest_genes_passing)
}

summarise_gene_expression <- function(
  object,
  genes,
  celltypes,
  condition_values
) {
  genes <- intersect(unique(genes), rownames(object))

  if (length(genes) == 0) {
    return(tibble())
  }

  expression_matrix <- GetAssayData(
    object,
    assay = "RNA",
    layer = "data"
  )

  purrr::map_dfr(celltypes, function(celltype) {
    purrr::map_dfr(condition_values, function(condition) {
      cells <- get_cells(
        object = object,
        celltype = celltype,
        condition = condition
      )

      if (length(cells) == 0) {
        return(tibble())
      }

      selected_expression <- expression_matrix[
        genes,
        cells,
        drop = FALSE
      ]

      tibble(
        celltype = celltype,
        condition = condition,
        gene = genes,
        n_cells = length(cells),
        n_positive = as.numeric(
          Matrix::rowSums(selected_expression > 0)
        ),
        pct_expressed = as.numeric(
          Matrix::rowMeans(selected_expression > 0)
        ),
        mean_expression = as.numeric(
          Matrix::rowMeans(selected_expression)
        )
      )
    })
  })
}

# STEP 1: Identify Cytokines Altered in the Senders
message("\n--- Step 1: Calculating Sender DEGs ---")

sender_de_list <- list()

for (sender in sender_celltypes) {
  message("Testing sender: ", sender)

  de <- tryCatch({
    FindMarkers(
      eae_obj,
      ident.1 = "EAE_ctrl",
      ident.2 = "EAE_ko",
      group.by = "eae_condition",
      subset.ident = sender,
      logfc.threshold = 0,
      min.pct = 0,
      assay = "RNA",
      verbose = FALSE
    ) %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::mutate(
        sender = sender,
        direction = ifelse(
          avg_log2FC > 0,
          "Higher in WT",
          "Higher in KO"
        ),
        condition_higher = ifelse(
          avg_log2FC > 0,
          "EAE_ctrl",
          "EAE_ko"
        )
      )
  }, error = function(e) {
    message("  -> Skipped due to error: ", e$message)
    tibble()
  })

  sender_de_list[[sender]] <- de
}

sender_de_all <- dplyr::bind_rows(sender_de_list)

altered_sender_chemokines <- sender_de_all %>%
  dplyr::filter(
    gene %in% cytokine_ligands,
    p_val_adj < padj_cutoff,
    abs(avg_log2FC) >= lfc_cutoff
  ) %>%
  dplyr::arrange(
    sender,
    p_val_adj,
    dplyr::desc(abs(avg_log2FC))
  )

write_csv(altered_sender_chemokines, file.path( output_dir, "01_altered_sender_chemokines.csv"))


# STEP 2: Map Altered Ligands to Expressed Receptors
message("\n--- Step 2: Mapping Receptors to Receivers ---")

target_receptors <- lr_network %>%
  dplyr::filter(
    from %in% altered_sender_chemokines$gene
  ) %>%
  dplyr::pull(to) %>%
  unique()

target_receptors <- union(
  target_receptors,
  interest_receptors
)

sender_ligand_expression <- summarise_gene_expression(
  object = eae_obj,
  genes = unique(altered_sender_chemokines$gene),
  celltypes = sender_celltypes,
  condition_values = conditions
) %>%
  dplyr::rename(
    sender = celltype,
    ligand = gene,
    sender_n_cells = n_cells,
    sender_n_positive = n_positive,
    sender_pct_expressed = pct_expressed,
    sender_mean_expression = mean_expression
  )

receptor_expression <- summarise_gene_expression(
  object = eae_obj,
  genes = target_receptors,
  celltypes = receiver_celltypes,
  condition_values = conditions
) %>%
  dplyr::rename(
    receiver = celltype,
    receptor = gene,
    receiver_n_cells = n_cells,
    receptor_n_positive = n_positive,
    receptor_pct_expressed = pct_expressed,
    receptor_mean_expression = mean_expression
  )

plausible_paths <- altered_sender_chemokines %>%
  dplyr::select(
    sender,
    ligand = gene,
    ligand_direction = direction,
    condition_higher,
    ligand_log2FC = avg_log2FC,
    ligand_p_value = p_val,
    ligand_padj = p_val_adj
  ) %>%
  dplyr::inner_join(
    sender_ligand_expression,
    by = c(
      "sender",
      "ligand",
      "condition_higher" = "condition"
    ),
    relationship = "many-to-many"
  ) %>%
  dplyr::filter(
    ifelse(
      ligand %in% interest_ligands,
      sender_pct_expressed >= interest_min_pct_expr,
      sender_pct_expressed >= min_pct_expr
    )
  ) %>%
  dplyr::inner_join(
    lr_network,
    by = c("ligand" = "from"),
    relationship = "many-to-many"
  ) %>%
  dplyr::rename(receptor = to) %>%
  dplyr::inner_join(
    receptor_expression,
    by = c(
      "receptor",
      "condition_higher" = "condition"
    ),
    relationship = "many-to-many"
  ) %>%
  dplyr::filter(
    ifelse(
      receptor %in% interest_receptors,
      receptor_pct_expressed >= interest_min_pct_expr,
      receptor_pct_expressed >= min_pct_expr
    )
  ) %>%
  dplyr::left_join(
    receiver_classification,
    by = "receiver"
  ) %>%
  dplyr::arrange(
    sender,
    receiver,
    ligand_padj,
    ligand,
    receptor
  )

View(plausible_paths)

message(
  sprintf(
    "Identified %d plausible ligand-receptor paths.",
    nrow(plausible_paths)
  )
)

View(plausible_paths)
write_csv(plausible_paths, file.path( output_dir, "02_plausible_ligand_receptor_paths.csv"))

# STEP 3: Calculate Receiver Differential Expression
message("\n--- Step 3: Calculating Receiver DEGs ---")

receiver_de_list <- list()

for (receiver in receiver_celltypes) {
  message("Testing receiver: ", receiver)

  receiver_de <- tryCatch({
    FindMarkers(
      eae_obj,
      ident.1 = "EAE_ctrl",
      ident.2 = "EAE_ko",
      group.by = "eae_condition",
      subset.ident = receiver,
      logfc.threshold = 0,
      min.pct = 0,
      assay = "RNA",
      verbose = FALSE
    ) %>%
      tibble::rownames_to_column("gene")
  }, error = function(e) {
    message("  -> Skipped due to error: ", e$message)
    tibble()
  })

  receiver_de_list[[receiver]] <- receiver_de
}

receiver_de_all <- purrr::imap_dfr(
  receiver_de_list,
  function(receiver_de, receiver) {
    if (nrow(receiver_de) == 0) {
      return(tibble())
    }

    receiver_de %>%
      dplyr::mutate(receiver = receiver) %>%
      dplyr::relocate(receiver)
  }
)

# STEP 4: Run Sender-Focused NicheNet
message("\n--- Step 4: Running Sender-Focused NicheNet ---")

comparison_directions <- tibble(
  program = c("Higher in WT", "Higher in KO"),
  condition_oi = c("EAE_ctrl", "EAE_ko"),
  condition_reference = c("EAE_ko", "EAE_ctrl")
)

nichenet_results <- list()
analysis_log <- list()
result_index <- 1L

for (sender in sender_celltypes) {
  for (receiver in receiver_celltypes) {
    receiver_de <- receiver_de_list[[receiver]]

    if (is.null(receiver_de) || nrow(receiver_de) == 0) {
      next
    }

    receiver_analysis_class <- receiver_classification$analysis_class[
      match(
        receiver,
        receiver_classification$receiver
      )
    ]

    for (direction_index in seq_len(nrow(comparison_directions))) {
      program <- comparison_directions$program[[direction_index]]
      condition_oi <- comparison_directions$condition_oi[[direction_index]]
      condition_reference <- comparison_directions$condition_reference[[direction_index]]

      message(
        "\n========================================\n",
        "Sender: ", sender,
        "\nReceiver: ", receiver,
        "\nProgramme: ", program,
        "\n========================================"
      )

      if (program == "Higher in WT") {
        geneset_oi <- receiver_de %>%
          dplyr::filter(
            p_val_adj < padj_cutoff,
            avg_log2FC >= lfc_cutoff
          ) %>%
          dplyr::pull(gene)
      } else {
        geneset_oi <- receiver_de %>%
          dplyr::filter(
            p_val_adj < padj_cutoff,
            avg_log2FC <= -lfc_cutoff
          ) %>%
          dplyr::pull(gene)
      }

      background_genes <- get_expressed_genes(
        object = eae_obj,
        celltype = receiver,
        condition = NULL,
        expression_threshold = min_pct_expr
      )

      background_genes <- intersect(
        background_genes,
        nichenet_target_universe
      )

      geneset_oi <- intersect(
        geneset_oi,
        background_genes
      )

      sender_expressed_genes <- get_eligible_sender_genes(
        object = eae_obj,
        celltype = sender,
        condition = condition_oi
      )

      receiver_expressed_genes <- get_eligible_receiver_genes(
        object = eae_obj,
        celltype = receiver,
        condition = condition_oi
      )

      potential_ligands <- lr_network %>%
        dplyr::filter(
          from %in% sender_expressed_genes,
          to %in% receiver_expressed_genes
        ) %>%
        dplyr::pull(from) %>%
        unique() %>%
        intersect(nichenet_ligand_universe)

      if (length(geneset_oi) < minimum_geneset_size) {
        message(
          "Skipped: fewer than ",
          minimum_geneset_size,
          " receiver genes."
        )

        analysis_log[[result_index]] <- tibble(
          sender = sender,
          receiver = receiver,
          analysis_class = receiver_analysis_class,
          program = program,
          condition_oi = condition_oi,
          condition_reference = condition_reference,
          n_geneset = length(geneset_oi),
          n_background = length(background_genes),
          n_potential_ligands = length(potential_ligands),
          status = "Skipped: insufficient receiver genes"
        )

        result_index <- result_index + 1L
        next
      }

      if (length(potential_ligands) == 0) {
        message("Skipped: no eligible sender ligands.")

        analysis_log[[result_index]] <- tibble(
          sender = sender,
          receiver = receiver,
          analysis_class = receiver_analysis_class,
          program = program,
          condition_oi = condition_oi,
          condition_reference = condition_reference,
          n_geneset = length(geneset_oi),
          n_background = length(background_genes),
          n_potential_ligands = 0L,
          status = "Skipped: no eligible sender ligands"
        )

        result_index <- result_index + 1L
        next
      }

      ligand_activities <- tryCatch({
        predict_ligand_activities(
          geneset = geneset_oi,
          background_expressed_genes = background_genes,
          ligand_target_matrix = ligand_target_matrix,
          potential_ligands = potential_ligands
        ) %>%
          dplyr::arrange(
            dplyr::desc(aupr_corrected)
          ) %>%
          dplyr::mutate(
            nn_rank = dplyr::row_number(),
            sender = sender,
            receiver = receiver,
            analysis_class = receiver_analysis_class,
            program = program,
            condition_oi = condition_oi,
            condition_reference = condition_reference,
            n_geneset = length(geneset_oi),
            n_background = length(background_genes),
            n_potential_ligands = length(potential_ligands)
          ) %>%
          dplyr::relocate(
            sender,
            receiver,
            analysis_class,
            program,
            condition_oi,
            condition_reference,
            n_geneset,
            n_background,
            n_potential_ligands
          )
      }, error = function(e) {
        message("Skipped due to NicheNet error: ", e$message)
        tibble()
      })

      if (nrow(ligand_activities) > 0) {
        result_key <- paste(
          sender,
          receiver,
          condition_oi,
          sep = "__"
        )

        nichenet_results[[result_key]] <- ligand_activities

        analysis_log[[result_index]] <- tibble(
          sender = sender,
          receiver = receiver,
          analysis_class = receiver_analysis_class,
          program = program,
          condition_oi = condition_oi,
          condition_reference = condition_reference,
          n_geneset = length(geneset_oi),
          n_background = length(background_genes),
          n_potential_ligands = length(potential_ligands),
          status = "Completed"
        )
      } else {
        analysis_log[[result_index]] <- tibble(
          sender = sender,
          receiver = receiver,
          analysis_class = receiver_analysis_class,
          program = program,
          condition_oi = condition_oi,
          condition_reference = condition_reference,
          n_geneset = length(geneset_oi),
          n_background = length(background_genes),
          n_potential_ligands = length(potential_ligands),
          status = "Skipped: NicheNet returned no results"
        )
      }

      result_index <- result_index + 1L
    }
  }
}

all_nichenet_activities <- dplyr::bind_rows(
  nichenet_results
)

analysis_log <- dplyr::bind_rows(
  analysis_log
)

write_csv(all_nichenet_activities, file.path(output_dir, "04_all_nichenet_activities.csv"))

write_csv(analysis_log, file.path(output_dir, "04_nichenet_analysis_log.csv"))

# STEP 5: Prioritise Interactions Supported at All Levels
message("\n--- Step 5: Prioritising Mechanistic Hypotheses ---")

if (nrow(all_nichenet_activities) > 0) {
  final_interactions <- plausible_paths %>%
    dplyr::inner_join(
      all_nichenet_activities,
      by = c(
        "sender",
        "receiver",
        "ligand" = "test_ligand",
        "ligand_direction" = "program",
        "condition_higher" = "condition_oi",
        "analysis_class"
      ),
      relationship = "many-to-many"
    ) %>%
    dplyr::arrange(
      analysis_class,
      sender,
      receiver,
      nn_rank,
      ligand_padj
    ) %>%
    dplyr::select(
      analysis_class,
      sender,
      receiver,
      direction = ligand_direction,
      condition = condition_higher,
      ligand,
      ligand_log2FC,
      ligand_padj,
      sender_pct_expressed,
      sender_mean_expression,
      receptor,
      receptor_pct_expressed,
      receptor_mean_expression,
      nn_auroc = auroc,
      nn_aupr = aupr,
      nn_aupr_corrected = aupr_corrected,
      nn_pearson = pearson,
      nn_rank,
      n_geneset,
      n_background,
      n_potential_ligands
    )
} else {
  final_interactions <- tibble()

  message(
    "Warning: no NicheNet ligand activities ",
    "were successfully calculated."
  )
}


# Save Main Outputs

write_csv(
  sender_de_all,
  file.path(
    output_dir,
    "01_all_sender_differential_expression.csv"
  )
)

write_csv(
  altered_sender_chemokines,
  file.path(
    output_dir,
    "01_sender_altered_chemokines.csv"
  )
)

write_csv(
  sender_ligand_expression,
  file.path(
    output_dir,
    "01_sender_ligand_expression.csv"
  )
)

write_csv(
  receptor_expression,
  file.path(
    output_dir,
    "02_receiver_receptor_expression.csv"
  )
)

write_csv(
  receiver_de_all,
  file.path(
    output_dir,
    "03_receiver_differential_expression.csv"
  )
)

write_csv(
  final_interactions,
  file.path(
    output_dir,
    "05_mechanistic_prioritised_interactions.csv"
  )
)

saveRDS(
  list(
    parameters = list(
      conditions = conditions,
      sender_celltypes = sender_celltypes,
      receiver_celltypes = receiver_celltypes,
      primary_receiver_celltypes = primary_receiver_celltypes,
      exploratory_receiver_celltypes = exploratory_receiver_celltypes,
      min_pct_expr = min_pct_expr,
      interest_min_pct_expr = interest_min_pct_expr,
      interest_ligands = interest_ligands,
      interest_receptors = interest_receptors,
      lfc_cutoff = lfc_cutoff,
      padj_cutoff = padj_cutoff
    ),
    cell_counts = cell_counts,
    sender_de = sender_de_all,
    altered_sender_chemokines = altered_sender_chemokines,
    sender_ligand_expression = sender_ligand_expression,
    receptor_expression = receptor_expression,
    plausible_paths = plausible_paths,
    receiver_de = receiver_de_all,
    all_nichenet_activities = all_nichenet_activities,
    analysis_log = analysis_log,
    final_interactions = final_interactions
  ),
  file.path(
    output_dir,
    "mechanistic_nichenet_full_results.rds"
  )
)