# Krzak, Willis et al 2025

## Succinate-SUCNR1 Signaling in Microglia Fuels Chronic CNS Inflammation

Grzegorz Krzak<sup>†</sup>, Cory M. Willis<sup>†</sup>, Evridiki Asimakidou<sup>†</sup>, Rafael Kollyfas, Maria Repollés-de-Dalmau, Nathan Whitty, Liviu Pirvan, Cristian Bulgaru, Julie A. Reisz, Rosana-Bristena Ionescu, Gregory Jordan, Monica Emili Garcia-Segura, Alexandra M. Nicaise, Samira Enache, Arianna Ghia, Ivan Lombardi, Berfin Barlas, David Rowitch, Gabriel Balmus, Angelo D’Alessandro, Sonia Fernández-Veledo, Irina Mohorianu, Stefano Pluchino, and Luca Peruzzotti-Jametti

† These authors contributed equally.

Myeloid cells play a pivotal role in persistent central nervous system (CNS) inflammation. Succinate, a tricarboxylic acid cycle metabolite, is a major regulator of myeloid cell functions through intracellular and extracellular signaling pathways. In this study, we aimed to unravel how extracellular succinate signaling via its succinate receptor 1 (SUCNR1) modulates myeloid cell function during CNS inflammation. We revealed that SUCNR1 exerts differential effects on myeloid cells: it blunts the pro-inflammatory polarization of CNS-infiltrating macrophages while it enhances that of CNS-resident microglia. Accordingly, using microglia-specific SUCNR1 knockout mice, we obtained a significant reduction of pro-inflammatory cytokines signaling and cerebrospinal fluid succinate in a multiple sclerosis (MS)-like disease mouse model, which was coupled with decreased immune infiltration and axonal damage. We further validated the pro-inflammatory effects of microglial SUCNR1 in post-mortem MS brains and in human microglia in vitro, thus establishing its critical role in regulating myeloid responses during CNS inflammation.

## Repository structure

- `SC_mouse/` – methods scripts for the ex vivo mouse spinal-cord scRNA-seq analyses (clustering annotation, composition statistics, differential expression, enrichment, and targeted displays). See [`SC_mouse/README.md`](SC_mouse/README.md) for install and run order. Large Seurat / ClustAssess RDS inputs are not bundled; place them under `SC_mouse/data/` as described there.
- `Scripts/` – R scripts included in the original submission for cluster-stability assessment with ClustAssess, BulkAnalyseR analysis, and ShinyCell visualization.
- `CYTOF/` – mass-cytometry workflow for preprocessing and balanced downsampling, UMAP and FlowSOM clustering, cell-type annotation, marker-expression summaries, and differential metacluster-abundance testing.
  - `CYTOF/CyTOF_data/` – input FCS files for the global *Sucnr1* knockout and the acute and chronic microglia-specific *Sucnr1* knockout datasets.
- `Cellchat_nichenet/` – preparation and analysis of the ex vivo scRNA-seq dataset with CellChat and NicheNet to compare WT and knockout cell–cell communication and prioritize microglial ligand–receptor–receiver hypotheses.
- `Human_regulons/` – analysis of SUCNR1 expression and transcription-factor regulon activity in the Absinta, Schirmer, and MacNair human single-cell or single-nucleus RNA-seq datasets, including pySCENIC input preparation and downstream visualization.
  - `Human_regulons/scenic_scripts/` – shell scripts for setting up pySCENIC, downloading human reference resources, and running the three dataset-specific pySCENIC workflows.
- `Spatial/` – spatial transcriptomics preprocessing and neighborhood-based analysis of a SUCNR1-associated gene panel across tissue niches, including expression summaries and statistical comparisons.
