library(ggplot2)
library(patchwork)
library(scales)
library(tidyverse)

output_dir <- "/iss-scratch/CoreBioinformatics/rk720/nichenet_stuff/mechanistic_nichenet_corrected"


# Load results
altered_sender_chemokines <- read_csv(
  file.path(output_dir, "01_sender_altered_chemokines.csv"),
  show_col_types = FALSE
)

final_interactions <- read_csv(
  file.path(output_dir, "05_mechanistic_prioritised_interactions.csv"),
  show_col_types = FALSE
)

all_nichenet_activities <- read_csv(
  file.path(output_dir, "04_all_nichenet_activities.csv"),
  show_col_types = FALSE
)

receptor_expression <- read_csv(
  file.path(output_dir, "02_receiver_receptor_expression.csv"),
  show_col_types = FALSE
)

sender_ligand_expression <- read_csv(
  file.path(output_dir, "01_sender_ligand_expression.csv"),
  show_col_types = FALSE
)

network_dir <- "/iss-scratch/CoreBioinformatics/rk720/nichenet_stuff/useful_rds"

receiver_de <- read_csv(
  file.path(output_dir, "03_receiver_differential_expression.csv"),
  show_col_types = FALSE
)

ligand_target_matrix <- readRDS(
  file.path(network_dir, "ligand_target_matrix.rds")
)

top_n_activity_ligands <- 10
top_n_target_ligands <- 5
top_n_target_genes <- 20


# Plot 1: Altered cytokines in microglial senders
plot_sender_de <- altered_sender_chemokines %>%
  dplyr::mutate(
    neg_log10_padj = -log10(p_val_adj),
    direction = factor(
      direction,
      levels = c("Higher in WT", "Higher in KO")
    )
  ) %>%
  ggplot(
    aes(
      x = avg_log2FC,
      y = reorder(gene, avg_log2FC),
      size = neg_log10_padj,
      colour = direction
    )
  ) +
  geom_point(alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(
    ~ sender,
    scales = "free_y"
  ) +
  scale_size_continuous(
    range = c(2.5, 8)
  ) +
  labs(
    title = "Cytokines and chemokines altered in microglial senders",
    subtitle = "Positive log2FC indicates higher expression in EAE control",
    x = "WT versus KO log2 fold-change",
    y = NULL,
    size = "-log10 adjusted p-value",
    colour = NULL
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )


# Prioritised interactions with WT and KO columns
interaction_plot_data2 <- final_interactions %>%
  dplyr::filter(
    analysis_class == "Primary",
    nn_aupr_corrected > 0,
    nn_rank <= 20 | ligand == "Cxcl2"
  ) %>%
  dplyr::select(
    sender,
    receiver,
    direction,
    ligand,
    receptor,
    nn_aupr_corrected,
    nn_rank
  ) %>%
  dplyr::distinct() %>%
  dplyr::left_join(
    receptor_expression %>%
      dplyr::select(
        receiver,
        receptor,
        expression_condition = condition,
        receptor_pct_expressed,
        receptor_mean_expression
      ),
    by = c("receiver", "receptor")
  ) %>%
  dplyr::filter(!is.na(expression_condition)) %>%
  dplyr::mutate(
    sender = stringr::str_replace_all(sender, "_", " "),
    receiver = stringr::str_replace_all(receiver, "_", " "),
    expression_condition = dplyr::recode(
      expression_condition,
      EAE_ctrl = "WT",
      EAE_ko = "KO"
    ),
    expression_condition = factor(
      expression_condition,
      levels = c("WT", "KO")
    ),
    interaction = paste0(ligand, " → ", receptor),
    interaction = forcats::fct_reorder(
      interaction,
      nn_aupr_corrected
    )
  )

plot_lr2 <- interaction_plot_data2 %>%
  ggplot(
    aes(
      x = expression_condition,
      y = interaction,
      size = receptor_pct_expressed,
      colour = receptor_mean_expression
    )
  ) +
  geom_point(alpha = 0.9) +
  facet_wrap(
    ~ sender + receiver,
    scales = "free_y",
    ncol = 2
  ) +
  scale_size_continuous(
    labels = scales::percent_format(accuracy = 1),
    range = c(2.5, 9)
  ) +
  labs(
    title = "Prioritised cytokine–receptor interactions",
    subtitle = "Bubbles show receptor expression separately in WT and KO",
    x = NULL,
    y = "Ligand → receptor",
    size = "Receiver cells expressing receptor",
    colour = "Mean receptor expression"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "bottom"
  )


# Targeted CXCL2–CXCR2–ACKR1 expression check
cxcl2_expression <- bind_rows(
  sender_ligand_expression %>%
    dplyr::filter(
      ligand == "Cxcl2",
      sender %in% c(
        "Microglia_Homeostatic",
        "Microglia_DAM"
      )
    ) %>%
    dplyr::transmute(
      celltype = sender,
      condition,
      gene = ligand,
      pct_expressed = sender_pct_expressed,
      mean_expression = sender_mean_expression
    ),

  receptor_expression %>%
    dplyr::filter(
      receptor %in% c("Cxcr2", "Ackr1"),
      receiver %in% c(
        "Neutrophils",
        "Endothelial_Cells",
        "Monocyte_derived_cells",
        "Oligodendrocytes"
      )
    ) %>%
    dplyr::transmute(
      celltype = receiver,
      condition,
      gene = receptor,
      pct_expressed = receptor_pct_expressed,
      mean_expression = receptor_mean_expression
    )
) %>%
  dplyr::mutate(
    celltype = stringr::str_replace_all(celltype, "_", " "),
    condition = recode(
      condition,
      EAE_ctrl = "EAE control",
      EAE_ko = "EAE KO"
    )
  )

plot_cxcl2_axis <- cxcl2_expression %>%
  ggplot(
    aes(
      x = gene,
      y = celltype,
      size = pct_expressed,
      colour = mean_expression
    )
  ) +
  geom_point(alpha = 0.9) +
  facet_wrap(
    ~ condition
  ) +
  scale_size_continuous(
    labels = percent_format(accuracy = 1),
    range = c(1, 11)
  ) +
  labs(
    title = "Targeted CXCL2–CXCR2–ACKR1 axis",
    subtitle = "CXCL2 is enriched in WT DAM microglia; CXCR2 is detected in neutrophils and ACKR1 mainly in endothelial cells",
    x = NULL,
    y = NULL,
    size = "Cells expressing gene",
    colour = "Mean expression"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )


# Save plots
ggsave(
  file.path(output_dir, "plot_sender_cytokine_changes.pdf"),
  plot_sender_de,
  width = 8,
  height = 7
)


ggsave(
  file.path(output_dir, "plot_prioritised_ligand_receptor_interactions.pdf"),
  plot_lr2,
  width = 9,
  height = 9
)


ggsave(
  file.path(output_dir, "plot_cxcl2_targeted_axis.pdf"),
  plot_cxcl2_axis,
  width = 9,
  height = 6
)


# NicheNet ligand activities for final prioritised interactions
#
# Receptor rows are collapsed because NicheNet activity is
# calculated for the ligand-receiver programme, not separately
# for each receptor.

final_activity_plot_data <- final_interactions %>%
  dplyr::filter(analysis_class == "Primary") %>%
  dplyr::group_by(
    sender,
    receiver,
    direction,
    ligand
  ) %>%
  dplyr::summarise(
    nn_aupr_corrected = dplyr::first(nn_aupr_corrected),
    ligand_log2FC = dplyr::first(ligand_log2FC),
    receptors = paste(sort(unique(receptor)), collapse = ", "),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    sender = stringr::str_replace_all(sender, "_", " "),
    receiver = stringr::str_replace_all(receiver, "_", " "),
    direction = factor(
      direction,
      levels = c("Higher in WT", "Higher in KO")
    ),
    panel = paste(receiver, direction, sep = "\n"),
    ligand_panel = paste(panel, ligand, sep = "___"),
    ligand_panel = forcats::fct_reorder(
      ligand_panel,
      nn_aupr_corrected
    )
  )

write.csv(final_interactions[final_interactions$ligand == "Ccl2",], paste0(output_dir, "/test.csv"))

plot_final_nichenet_activity <- final_activity_plot_data %>%
  ggplot(
    aes(
      x = nn_aupr_corrected,
      y = ligand_panel,
      fill = direction
    )
  ) +
  geom_col(width = 0.7) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ panel,
    scales = "free_y",
    ncol = 2
  ) +
  scale_y_discrete(
    labels = function(x) sub("^.*___", "", x)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.1, 0.25))
  ) +
  labs(
    title = "NicheNet activities of mechanistically prioritised ligands",
    subtitle = "Ligands passed sender differential-expression and receptor-expression filtering",
    x = "Corrected AUPR",
    y = NULL,
    fill = NULL
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    output_dir,
    "plot_final_interaction_nichenet_activities.pdf"
  ),
  plot_final_nichenet_activity,
  width = 10,
  height = 10
)