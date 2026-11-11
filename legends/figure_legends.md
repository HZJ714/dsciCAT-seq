Figure legends

Figures are submitted as separate files, one per figure (PDF or TIFF, 16.0-18.0 cm wide; Figure S6, 17.0 x 18.9 cm; Figure S7, 17.0 x 14.28 cm). The Figure S6 and Figure S7 legends below match the final seven-panel (A–G) layouts. The Graphical Abstract is supplied separately (20.5 x 11.6 cm).

Figure 1 | dsciCAT-seq jointly profiles transcriptome and chromatin accessibility in the same nucleus, reaching 100,000-nucleus input with sensitivity comparable to established chemistries.

A, Experimental and computational workflow of dsciCAT-seq. The Tn5 and reverse-transcription reactions of one nucleus carry matched 10-nt pre-indices; the droplet barcode is added on the beads afterwards, and the UMI and sequencing index are separate from both.

B, Species mixing of dsciCAT-seq-100k with 16 indices (RNA top, ATAC bottom; dots coloured by species). Of 42,816 called barcodes, 19,397 are human-only, 19,597 mouse-only and 3,822 (8.93%) carry reads of both species, so the dual-modality singlet set is 38,994 (91.07%).

C, Genome-browser tracks over one locus per species (human MIR17HG in HEK293T, mouse Hes7 in NIH3T3) for dsciCAT-seq-50k, dsciCAT-seq-100k and HT-scCAT-seq; tile rows show the 120 cells with the most cuts in each window, on a y-scale shared within each species block.

D, Box plots of UMIs (top) and detected genes (bottom) per cell, log10 scale. Methods marked * (SUM-seq, Paired-seq, sci-CAR) use released count matrices and are not depth-matched; each method/cell-line combination is one technical library, so cross-method comparison is descriptive. Numbers of cells used per library and modality are listed in Table 1.

E, Violin plots of unique fragments per nucleus, same layout as D (white dot = median, black bar = IQR); the two modalities are filtered independently. ‡, the sci-CAR NIH3T3 ATAC column is quantified on a human-only peak set and reflects the reference, not sci-CAR performance; the RNA column of the same accession is normal. Numbers of nuclei used per library and modality, including the unfiltered sci-CAR matrix, are listed in Table 1.

F, Cross-method transcriptome similarity: Pearson correlation of pseudobulk profiles over genes detected in every method (log2 CPM + 1; colour scale 0.2–1.0). Row and column swatches repeat the method colours of D and E.

Figure 2 | dsciCAT-seq captures paired RNA and ATAC profiles across the adult mouse brain and, integrated with public 10x Multiome and HT-scCAT-seq data, resolves 33 cell types including oligodendrocyte-lineage subtypes.

A–C, UMAPs of the 33 cell types of 70,643 nuclei in RNA (A), ATAC (B) and WNN (C) space; numbers indicate clusters in legend order (colour table reused in Figures S2 and S3).

D, Dot plot of 48 marker genes across the 33 cell types; Dot size gives the fraction of expressing cells and colour the average expression level, scaled within each gene row.

E, Genome-browser tracks of pseudobulk ATAC coverage for Qk, Gpr17, Plp1, Slc32a1 and Slc17a7, aligned row by row with D (each row on its own maximum).

Figure 3 | Raising nucleus input from 50,000 to 100,000 recovers the rare COP/NFOL lineage population at sufficient cell number, while the lineage structure itself is unchanged.

A–C, WNN UMAPs for HT-scCAT-seq (A, 4,029 nuclei), dsciCAT-seq-50k (B, 17,093) and dsciCAT-seq-100k (C, 30,350); each panel is one library, not a replicate. Circled COP/NFOL clusters are expanded into Gpr17 feature plots (B, 60 nuclei; C, 133 nuclei; sub-clusters: Figure S3I, J).

D, UMAP of the six subtypes of the oligodendrocyte-lineage subset (5,977 nuclei; RNA-only embedding, distinct from A–C).

E, Dot plot of 17 lineage markers across the six subtypes; size = percentage of expressing cells, colour = average expression.

F, Genome-browser tracks of peak coverage for the subtype markers Bmp4 (COP), Synpr (MFOL), Anln (MOL2), Nrxn3 (MOL3) and Grm3 (MOL6), aligned one-to-one with the rows of E (each row on its own maximum).

Figure 4 | Chromatin–RNA dynamics along the oligodendrocyte lineage.

(A) scMEGA joint RNA–ATAC uniform manifold approximation and projection (UMAP) of 5,977 oligodendrocyte-lineage cells, colored by pseudotime (0–100): complex/oligodendrocyte precursor-like–newly formed oligodendrocyte-like (COP/NFOL) cells at one pole, myelinating-fiber oligodendrocyte-like (MFOL) in transit, and the mature subtypes MOL2, MOL3, MOL3/6 and MOL6 at the other. Panels A and B are drawn on the RNA-based embedding produced by scMEGA; the weighted-neighbourhood (WNN) UMAP shown in Fig 3 (six-state version: Fig 3D) was not used for Fig 4, and the cell order of the heat maps in (B) is the pseudotime order of (A), not an independent accessibility embedding.

(B) Paired heat maps of chromatin accessibility (top; 15,607 peak-to-gene rows) and gene expression (bottom) for the same 5,977 cells ordered along the pseudotime axis of (A) (vertical axis, 0 at the bottom to 100 at the top; the same pseudotime scale as the color bar of (A), shown here as plain ticks rather than a second color bar); colors are per-row z-scores clipped at ±2 (top scale bars). The horizontal axis is the 15,607 peak-to-gene links, ordered by the row dendrogram of the combined z-score matrix, so columns are features rather than successive time points; it is not a genomic coordinate. Fourteen stage markers are annotated above the accessibility panel with leader lines (early Tcf7l2, Bmp4 and Olig2; transitional Gpr17, Mobp, Opalin, Mog, Car2, Apod, Mal and Olig1; mature Plp1, Mag and Mbp; gene symbols in italics denote mouse orthologs).

(C) Summary of kinetic-model statistics for the 1,809 fitted genes (model assignment is shown in Fig S4C: 1,798 = 99.4% model 1, the chromatin-dependent three-phase model; 11 = 0.6% model 2, the two-phase model, including Plp1 and Tcf7l2). The horizontal axis of the left panel is the signed 2ΔlogL of the chromatin–RNA decoupling test — twice the log-likelihood difference between the fit with and without a decoupling phase — so positive values favour the decoupled model: 372 genes (20.6%) are decoupled at a false discovery rate [FDR] < 0.05, whereas 1,434 of 1,809 (79.3%) are negative because adding the phase worsens the fit and their one-sided P is then 1 by construction, as printed. 208 genes return an FDR of exactly 0 through double-precision underflow of the χ² tail and are better read as P < 1 × 10⁻³⁰⁰, so −log₁₀(P) saturates at 300 and the maximum 2ΔlogL is 3.3 × 10⁴. Induction timing: 93.6% of the 1,798 model-1 genes have a first switch time t_sw1 < 0, i.e. chromatin accessibility opens before the expression switch, here termed priming (descriptive temporal precedence, not causal evidence). Switch-interval lengths (as a proportion of the 0–20 gene-time window) for the 114 high-confidence velocity genes: primed n = 17, coupled-on n = 17, decoupled n = 18, coupled-off n = 10; each box is the subset of the 114 with a well-defined interval within the window (total n = 62). Decoupled intervals were longer than all other interval classes combined (median 0.92 vs 0.33 pooled across the other three classes, on the 0–20 gene-time window scale; Mann–Whitney U = 640, n = 18 versus 44 of the 62 valid intervals; one-sided P = 8.0 × 10⁻⁵, two-sided P = 1.6 × 10⁻⁴; the one-sided direction — decoupled longer — was specified before testing; bootstrap 95% CI of the median difference [0.31, 0.87]). The decoupled class in this interval panel is a switch-interval classification within the 114-gene set and is not the same set as the 372 genes with a positive signed 2ΔlogL; the two intersect in only 32 genes.

(D) MultiVelo dynamical fits for three representative genes: Anln (model 1, r² = 0.58, all three fitted switches inside the window, 114-gene high-confidence set), Cemip2 (model 1, r² = 0.83, in-window switches but outside the high-confidence set) and Tcf7l2 (model 2, r² = 0.90, transcription-led with flat fitted accessibility), showing accessibility, unspliced (Mu) and spliced (Ms) RNA against the fitted gene-time axis (0–20); dashed lines mark the three fitted switch times. Accessibility axes are displayed to the 99th percentile of the per-cell fitted values (≤ 1% of cells above the top edge) to keep the fitted curves legible. Genes whose fitted switches lie at or beyond the window edge (for example Tnr and Plp1) are shown descriptively in Figure S4D.

Figure 5 | SCENIC+ reconstructs subtype-specific eRegulons and identifies a Mef2c-centred regulatory program in the MOL3 subtype.

A, Heatmap with dots of the activity of 20 subtype-specific eRegulons across the six subtypes. Row labels give the motif-driving transcription factor and the number of regulatory regions; asterisks mark extended eRegulons.

B, Gene regulatory network of the MOL3 subtype. Edges are regulator-to-target links taken from the eRegulons (TF2G table only; no peak-region nodes are drawn), so direction is read from the node type: enlarged diamonds with bold labels are the regulators named by the motif-to-transcription-factor annotation of SCENIC+, which assigns the driving motif to a gene name (Mbnl2, an RNA-binding splicing regulator, is one such label), and smaller diamonds are target genes. Edge colour identifies the source eRegulon and edge width increases linearly with TF2G importance (0.61-4.24 in this panel; 1.0-4.0 pt on the source canvas, 0.36-1.44 pt as printed). Gene node colour gives the MOL3 log2 fold-change, blue for low and red for high. All 44 retained links are predicted activators; no repressor eRegulon passed the correlation filter. A single-row graphical key (regulator diamond, arrow, target diamond, and the thinnest and thickest sample edges) is printed along the top of the panel, level with the panel letter; every other encoding is defined in this legend only. Links are inferred from motif enrichment and accessibility-expression correlation, not from direct binding assays.

C, Scatter plots of eRegulon activity within one technology (y) against dsciCAT-seq six-subtype activity (x, from that comparison’s own SCENIC+ run) for the 145 eRegulons shared by both. The two arms are 10x Multiome (four batches) and HT-scCAT-seq (four libraries), both integrated with Harmony. Dashed lines indicate identity; rho and P are from Spearman correlation.

Figure 6 | Multiple sclerosis genetic-risk scores are highest in mature oligodendrocyte subtypes, and the mature-oligodendrocyte-versus-precursor contrast replicates between the expression and chromatin channels of the same nuclei.

(A) Single-cell disease-risk score (scDRS) distribution for the MS GWAS gene set across six oligodendrocyte-lineage subtypes (5,962 scored nuclei of 5,977; COP/NFOL n = 203, MFOL n = 490, MOL2 n = 1,016, MOL3 n = 688, MOL3/6 n = 991, MOL6 n = 2,574). Boxes show the median and interquartile range, whiskers the 5th–95th percentiles and dots individual nuclei; colours are the manuscript-wide ArchR stallion subtype palette. The exact Mann–Whitney U (MWU) P value against the reference group COP/NFOL is printed above each box (MOL6 6.9 × 10⁻¹⁰, MOL2 1.3 × 10⁻⁸, MOL3/6 1.3 × 10⁻⁶, MFOL 5.0 × 10⁻⁵, MOL3 0.039). A more conservative group-level permutation control (scdrs_group, 1,000 control gene sets) reached nominal P < 0.05 only for MOL3/6 (0.008) and MOL6 (0.019), and no subtype survived false-discovery-rate correction; the score is therefore interpreted as a ranking statistic rather than an absolute risk measure.

(B) The same normalized score (norm_score) displayed as the mean per UMAP bin (hexagonal binning, ≥5 nuclei per bin) on the weighted-neighbourhood (WNN) joint RNA–ATAC UMAP of Fig 3 (the UMAP coordinates of panels B, S6C and S6F are that embedding’s axes), on a light-blue–white–red diverging scale centered at 0 with a symmetric range of ±0.74 (1st–99th percentile of the bin means). The COP/NFOL precursor tail at the lower left is light blue, whereas the main oligodendrocyte cloud is red-tinted, with the highest values in the mature subtypes MOL2 and MOL6, as annotated. Inference is based on the cell-level scores in (A); the binned view is a display choice only. The unbinned per-nucleus view and the cross-chemistry replication of the same contrast are shown in Figure S6F.

(C) The 74 GWAS-derived gene sets (MAGMA 10 kb windows, 1,000 genes each): the mature-oligodendrocyte-versus-precursor contrast Δ1 (median MOL6 − median COP/NFOL) computed in the gene-expression channel (y-axis) against the same contrast in the gene-score (chromatin-accessibility) channel (x-axis) of the same 5,962 nuclei. Spearman ρ = 0.819 (permutation P = 1 × 10⁻⁴, add-one corrected; 0 of 10,000 null draws, seed 20261104, reached the observed ρ). A trait-level control — the correlation of the two channels’ absolute scores — is ρ = 0.029 (P = 0.81), and ρ = 0.679 after removing the trait absolute level, so the lineage contrast, not the absolute score level, is what replicates between channels. Named traits give their within-channel rank: MS (RNA #1, gene score #18), schizophrenia (RNA #55, gene score #38) and height (RNA #6, gene score #1). The two axes are scaled independently, so the dashed identity line is not at 45° on the page.

Figure 7 | The mature-oligodendrocyte-versus-precursor risk contrast replicates in independent human MS cohorts and in mouse EAE, and is specific to the oligodendrocyte-maturation axis.

(A) Unit-level Δ1 = median scDRS score of mature oligodendrocytes − median score of OPC/COP precursors, computed within each donor (GSE277430, GSE118257, GSE231585), within each library (GSE301908) and within each mouse library (GSE250589, spinal cord, FACS GFP+ from Sox10-cre animals, lesion phase as labelled). Dots are units, bars are group medians, colours are MS (red), EAE (orange) and control (dark blue). The number above each column is the Mann–Whitney rank-sum P for that column (MS or EAE versus control); the “P =” prefix is omitted in the panel for legibility at the printed column width. Values: GSE277430 0.500, GSE118257 1.000, GSE301908 0.091, GSE231585 0.393; mouse Early 0.100, Peak 0.333, Late 0.857.

(B) The same contrast ranked across 74 GWAS gene sets (MAGMA, 10 kb windows, 1,000 genes each) in the four human cohorts. Rows are the 74 gene sets ordered top-down by the mean of their four within-cohort ranks; columns are cohorts; colour is the rank of cohort-level Δ1 within that cohort (1 = strongest mature-oligodendrocyte contrast, 74 = weakest). Because Δ1 tracks the size of each gene set’s regulatory domain in every cohort (Spearman ρ of Δ1 with ln n_peak = −0.84 / −0.69 / −0.78 / −0.87; n_peak counted on public human snATAC GSE231584, GRCh38, 121,505 deduplicated peaks, TSS ± 100 kb), the ranks shown are those obtained after residualising Δ1 on ln n_peak, so that a gene set is not ranked ahead of another merely for covering more chromatin. On this basis MS ranks 1 / 1 / 3 / 1, Alzheimer’s disease 18 / 49 / 27 / 12 and schizophrenia 68 / 74 / 42 / 71 (GSE277430 / GSE118257 / GSE301908 / GSE231585). The unadjusted ranking (MS 2 / 4 / 7 / 2), the ±15% matched-window control (MS 2 / 3 / 5 / 2) and the ρ strip are reported in Table S7. Ranks are always computed within a cohort and never from the transformed values in (C).

(C) The same nuclei re-projected onto two further cell-identity axes: astrocyte−immune (median astrocyte − median immune) and neuron−non-neuron (median neuron − median of all non-neuronal populations: glia, immune and vascular), beside the oligodendrocyte-maturation axis of (A) and (B). Values are robust z of Δ, standardised within cohort × axis by median and MAD, so the four tiles share one y-axis; violins are kernel densities of the 74 gene-set values in that cohort × axis, evaluated only inside the visible window (truncated at −4.7 / +5.5), while ranks and interquartile boxes use all 74. Red dots are MS and dark-blue dots schizophrenia, each labelled #n with its rank on the axis it sits on. The single asterisk marks that axis, not that cohort: GSE231585 contains 62,911 nuclei of which 717 (1.1%) are annotated neurons, 4–133 per sample, so under the frozen rule that a group must contribute at least 10 nuclei to a per-unit median, two of its 15 samples drop out of the neuron−non-neuron axis (13 units there, 14 on the other two axes). Nothing is subsampled; 717 is the annotated neuron pool of the cohort. The other cohorts’ neuron pools are 7–100-fold larger (GSE277430 76,436; GSE118257 5,031; GSE301908 ≤ 12,000 per lineage by the stratified draw), so no other tile carries a mark.

(D) The 18 approved or in-trial MS drug targets shown out of 80 scored. Left block: within-cohort log2 odds ratio of detection in mature oligodendrocytes versus OPC; middle block: OL lineage versus all non-lineage cells; the left block covers all four cohorts and the middle block the three cohorts that carry non-lineage annotations. Right block: the percentage of nuclei in which the gene is detected in GSE301908, shown separately for each cell population — a within-population detection rate, not a composition. The red box marks the single row that passes the pre-set rule (mature OL exceeds every non-lineage population by ≥1.0 percentage point and by ≥1.5-fold): S1PR5. Resolved per donor rather than pooled over nuclei (Fig S7G), the rule is met in 22 of 27 evaluable donors across the three cohorts that carry non-lineage populations (8/11 GSE301908, 7/8 GSE118257, 7/8 GSE231585; a population enters a donor’s comparison only if that donor contributes ≥ 30 nuclei of it), and all five failures are thin-comparator donors (31–410 nuclei in the population that beat or nearly beat OL). GSE301908’s 11 donors are 8 patients and 3 controls, and diagnosis does not separate the margin (pooled P = 0.34; per-cohort P = 0.571 / 1.000 / 0.376; Table S7); GSE277430 cannot carry the rule because its released object is an OL-enriched subset without a non-lineage comparator. No other target passes the rule in any evaluable donor. The three drug-target columns (named in a product label / stated mechanism of action / used in MS through this target) describe the drug-target relation itself and are not localisation evidence.

The statistical unit is the donor (library where a cohort is library-resolved), never the nucleus; n per unit, cohort sources, platforms and annotation sources are given in Table S7 and Table S6. scDRS 1.0.2 (score_cell, n_ctrl = 200, fixed seeds per cohort).

Figure S1 | Library quality of dsciCAT-seq is consistent across barcode numbers and throughput scales.

A, Bar and box plots of RNA-channel quality control of the eight libraries: cell numbers, UMIs and detected genes for HEK293T (left) and NIH3T3 (right).

B, Violin plots of unique fragments per cell, same layout as A.

C, Heatmaps of pairwise RNA similarity between the eight libraries in the two cell lines.

D, Stacked bars of dual-omic species mixing of all eight libraries, same layout as Figure 1B.

Figure S2 | Joint integration of three library chemistries from three separate studies, comprising ten libraries and analysed as four groups, supports cell-type composition comparison and the choice of sequencing depth.

A, UMAP of the 70,643 nuclei in WNN space, faceted by library.

B, Cell-type composition per library, in the same order as A; the two dsciCAT-seq libraries come from different animals and are not replicates.

C, Clustering consistency (NMI, ARI) against nuclei analysed; lines = mean over 20 downsampling replicates, bars = ±1 SD (computational, not biological, replicates).

D, Recall of the 37 cell types of the dsciCAT-seq-only object (59,962 nuclei; finer than the 33-type annotation) across eight throughput levels, ordered by abundance (rarest first). Red labels mark the 16 types with recall gain ≥ 0.5 at 40k versus 3k (paired test over 20 replicates, BH FDR < 0.05).

Figure S3 | Detection of COP/NFOL depends on the number of lineage nuclei analysed, and a pre-specified sub-clustering test bounds how far that structure can be claimed.

A, RNA UMAPs of HT-scCAT-seq (4,029), dsciCAT-seq-50k (17,093) and dsciCAT-seq-100k (30,350), as in Figure 3A–C.

B, ATAC UMAPs of the same libraries and cells.

C–E, Pearson correlation between average expression profiles of every cluster pair across libraries (5,597 shared variable genes); colour z-scored within each query column, printed r values give absolute strength. Red labels mark types the reference library does not resolve; the black box marks the COP/NFOL column and its best match. “300/300 bootstraps” = gene resamples in which the best match is unchanged (r = 0.813/0.803/0.986). Because both sides share markers and annotation framework, this is a consistency check, not independent validation.

F, COP/NFOL recovery against the number of lineage nuclei sub-clustered. Filled circles = full library, hatched = downsampled (fractions = successful draws out of 20), open crossed circle = no recovery. Candidates were called COP/NFOL by a rule fixed in advance (n ≥ 10, Bmp4 ≥ 0.25, Gpr17 ≥ 0.5, ≥30% Gpr17+ cells, C1ql1 ≤ 0.2, Ptprz1 ≤ 1.5). Cluster sizes refer to this sub-clustering and differ from Figure 3 (60/133) and the integrated lineage object (203).

G, Seven lineage markers for the candidate, OPC and mature clusters of the four conditions in F; red frames = clusters meeting the rule, the grey frame marks the n = 15 candidate failing the C1ql1 criterion.

H, COP/NFOL nuclei of both dsciCAT-seq libraries re-clustered with identical parameters (adjusted Rand index 1.00 against Figure 3C): the 50k library (60 nuclei) splits only at resolution 0.8 (32/28); the 100k library (133 nuclei) gives sub-0 (n = 81, 2,682 UMI) and sub-1 (n = 52, 4,292 UMI).

I, Frozen marker panel that took no part in defining COP/NFOL. Cdk1, Mki67, Top2a and Dcn (below the dashed line) are undetectable in all groups, excluding cycling or stromal contamination; sub-1 reaches the oligodendrocyte reference for myelin genes (Mobp, Mog, Cnp), sub-0 does not.

J, Recovery of the two-cluster split across a perturbation grid fixed before the data were seen (dark bars = whole-library runs; pale bars = 100 runs at 80% downsampling; green diamonds = adjusted Rand index against Figure 3C). Of the three pre-specified criteria, S1 is met (260 genes at BH FDR < 0.05), S2 is not met over the whole grid (41.9%, ARI 0.62) though met at resolution 0.3–0.5, and S3 is not met for nuclear RNA capture (2,682 versus 4,292 UMI; ATAC metrics clean). The split is therefore reported descriptively, not as two substates.

Figure S4 | Time axes, peak-to-gene links and gene-set structure behind Fig 4.

(A) MultiVelo latent time on the same UMAP embedding as Fig 4A (left; 5,977 cells, rank-coloured 0–1) and the cell-type composition of the latent-time axis (right; stacked subtype proportions in equal-width latent-time bins). Latent time is only weakly correlated with the scMEGA pseudotime of Fig 4A at the cell level (Spearman ρ = 0.20) and does not order the mature subtypes (MOL3/6 0.747, MOL2 0.665, MOL6 0.659), so the two time axes are not interchangeable and neither supports a monotone maturation gradient within the mature block.

(B) Properties of the 15,607 peak-to-gene links that underlie the Fig 4B accessibility heat map (2,486 genes, 13,010 peaks; scMEGA with correlation r ≥ 0.40, distance ≤ 250 kb and FDR ≤ 1 × 10⁻⁴): histograms of link correlation (median r = 0.56), peak–TSS distance (median 107 kb) and peaks per gene (median 5, maximum 27). All 331 S1 marker genes appear in the link table (orange). The axes are truncated at the link-selection cut-offs, so no enrichment, pass-rate or decay-rate statement is made from these histograms.

(C) Gene-set structure and programmes: the 1,809 fitted genes give two partially overlapping subsets—372 genes with a positive signed 2ΔlogL at FDR < 0.05 (Fig 4C) and 114 high-confidence velocity genes (82 coupled, 31 priming, 1 lagging shutdown)—that intersect in only 32 genes and are therefore not nested (funnel, top). Model-1 (1,798 genes, rows in k-means programme order) and model-2 (11 genes) smoothed spliced-RNA heat maps use the same latent-time cell order and per-gene z-scores; the four k-means programmes are C0 (n = 576, early: Vcan, Tnr, Gria2, Dnm3), C1 (n = 491, MFOL onset: Mbp, Ugt8a, St18, Rtn4), C2 (n = 433, mature sustained: Ptgds, Mal, Mobp, Cnp, Opalin) and C3 (n = 298, mid-transient: Apod, Trf, Anln, Kif1b). Direction examples (LRT here denotes the signed 2ΔlogL of Fig 4C): Nfasc (the single lagging-shutdown gene; model 2, LRT 700, FDR 2 × 10⁻¹⁵³) and Ndrg1 (coupled; model 1, LRT −6,755, FDR 1); dashed lines mark t_c50 and t_s50. Nfasc shows that a late accessibility switch relative to RNA (t_s50 0.4 < t_c50 9.8) does not by itself indicate RNA-first timing—the lagging call comes from the shutdown phase, not from the switch times.

(D) Three-dimensional spliced–unspliced–chromatin (s–u–c) phase portraits for eight genes spanning the categories above (Anln, Mbp, Mal and Car2, priming; Tcf7l2, coupled; Tnr, Plp1 and Ptgds, not in the 114-gene set), in raw model units; the black line is the fitted trajectory and magenta marks the switch times. In the panel headers, values labelled LRT are the signed 2ΔlogL statistic of Figure 4C (positive favours the decoupled model); FDR is reported as 1 where that statistic is negative and as <1e-300 where the χ² tail underflows. Each header gives the selected model, LRT, FDR and the switch times t_c50/t_s50; the view is identical in all eight panels.

Figure S5 | eRegulon selection, network structure, locus-level validation and down-sampling robustness behind Figure 5.

A, eRegulon specificity (eRSS) of the six subtypes on one shared axis; dashed line = eRSS 0.5. Red points mark each subtype’s top four eRegulons and everything above 0.5; full top-10 ranks in Table S5.

B, Gene regulatory networks of the MOL2, MOL3/6 and MOL6 subtypes, drawn as in Figure 5B (key lower left); all retained links are predicted activators.

C, ATAC insertion at the Nrxn3 locus per subtype (E1–E4 = Mef2c-eRegulon enhancers), with Nrxn3 expression per subtype (right); one-sided Mann-Whitney U, MOL3 versus the other five pooled: P < 1×10⁻⁷⁹.

D, eRegulon number, target-set Jaccard similarity and top-10 TF recovery against the number of cells analysed.

Figure S6 | Positional, gene-level and cross-chemistry evaluation of the oligodendrocyte-lineage MS risk signal.

(A) Four example loci (PHLDB1, PHGDH, PLLP and BCAS1), each spanning ±500 kb around the lead SNP of its MS genome-wide significant signal (GWAS Catalog, PMID 31604244, P ≤ 5 × 10⁻⁸, hg38): the black vertical line is the index SNP, the thick red segment the transcription start site (TSS) of the platform-discovered gene and grey segments other genes in the window; each panel states the SNP P value and the signed TSS–SNP distance, with the hit gene named below its TSS. The full list of the 38 platform genes whose TSS falls within a ±500-kb window, with per-gene distances, is in Table S6.

(B) Permutation test for the overlap rate. The observed 38 of 331 platform genes (11.5%) with a TSS inside a ±500-kb window is compared with a protein-coding background of 10.34% (2,071 of 20,020 GENCODE v45 basic genes). The null is recomputed from source (10,000 draws, seed 20261104, add-one correction): permutation P = 0.277; Fisher exact P = 0.467 (odds ratio 1.13) — no enrichment. The value of the overlap lies in the identity of the hit genes, not in the enrichment statistic; all 38 overlapping platform genes are protein-coding, so observed and permuted pools are of the same gene class.

(C) Bin-size robustness of the binned view of Fig 6B: the same nuclei, the same minimum of 5 nuclei per bin and the same ±0.74 colour scale, drawn as single cells and at hexbin grid sizes 15, 25 and 40 (147, 362 and 509 filled bins). The light-blue COP/NFOL tail and the red-tinted mature cloud are present at every resolution, so the pattern in Fig 6B is not an artifact of smoothing.

(D) Gene-level landing of the MS MAGMA top-1000 — the gene set actually used for scoring in Fig 6A/6B — across the six subtypes in both channels. Rows are the 223 top-1000 genes measurable on this platform (875 mapped minus 652 near-zero); the displayed quantity is the per-gene mature-block contrast Δ (the weighted MOL3/6-plus-MOL6 offset from the COP/NFOL reference), row-standardized within each channel. The set-level test against detection-matched non-MS background gene sets (decile-stratified permutation, 10,000 draws) is null: P = 0.096 (RNA) and P = 0.745 (gene score). The panel is deliberately unlabelled: a null panel with gene names on it would invite the reader to hunt for a pattern the test says is not there. Cells beyond the symmetric colour limit carry a small white triangle at the clipped edge of the saturated cell (shared legend, right).

(E) The cis-window set lands on MOL6. The 38 genes whose platform gene lies inside a lead-SNP ±500-kb window, same display and contrast as (D), give a positive set-level result in both channels: mean Δ = +0.145 (RNA) and +0.097 (gene score) against detection-breadth- and signal-matched backgrounds, z = +2.87 and +3.07, P = 0.0017 and 0.0007; the myelin/oligodendrocyte anchor gene set gives 0.0107 and 0.0080, and the 293 platform genes outside any window give 0.0127 and 0.0066 (leave-one-gene-out worst cases: 0.0127 dropping PRR5L, 0.0038 dropping KCTD13). Panels D and E are not nested: only 8 of the 38 window genes appear in the MS top-1000 at all (6 measurable here), so the two panels share 6 rows — two gene-level definitions of the same GWAS disagree, which is the finding these panels are there to show. The four named genes are the loci of panel A. BCAS1 is the one a reader may question: it is the progenitor-end gene of the four (top subtype COP/NFOL, Δ_MOL6 = −1.71, ranked 38 of 38; −1.14 in the independent human mature-vs-precursor expression contrast) and is not in the MS top-1000, so its negative Δ is not a counter-example to Fig 6A/6B.

(F) The same OPC-versus-mature-oligodendrocyte contrast that Fig 6A/6B states for the own cohort, put per nucleus into the three library chemistries and both channels: two rows are the channels, three columns the chemistries; the two-sided Mann–Whitney U (MWU) P value printed per panel is taken from Table S6. Own dsciCAT-seq: RNA P = 1.7 × 10⁻³³ (n = 1,488 OPC / 5,813 oligodendrocyte) and gene score P = 0.0022 (1,490 / 5,827); 10x Multiome: RNA P = 1.8 × 10⁻¹⁸ and gene score P = 3.1 × 10⁻⁴ (298 / 540); HT-scCAT-seq: RNA P = 6.8 × 10⁻⁶ and gene score P = 0.75 (220 / 293) — the HT chromatin channel is the one non-significant panel.

(G) The 74-trait expression/chromatin Δ1 pair across the three chemistries (the panel formerly shown as Fig 6D), re-cut for the supplement. Rows are the 74 GWAS gene sets ordered by ascending mean expression-channel size-adjusted rank (only the expression channel orders; the chromatin channel never enters the ordering). For each chemistry a Δ1 column (within-column z score, coolwarm, symmetric clip — absolute Δ1 values are never compared across chemistries) is shown beside a rank column (ranks are the only quantity compared across chemistries); only MS, schizophrenia and height are labelled. Size-adjusted cross-channel ρ = 0.896 (own), 0.824 (10x) and 0.630 (HT), each at the 10,000-draw permutation floor P = 1 × 10⁻⁴; bootstrap 95% confidence intervals 0.806–0.944, 0.700–0.902 and 0.451–0.768. The MS rank moves expression → chromatin: 7 → 37 (own), 4 → 49 (10x) and 11 → 41 (HT). All three chemistries meet the pre-specified ladder threshold (ρ_adj ≥ 0.60), so the size-adjusted cross-channel ladder reproduces in all three, never as a statement that “the ladder holds” beyond this sample.

Abbreviations: COP/NFOL, complex/oligodendrocyte precursor-like–newly formed oligodendrocyte-like; MOL2/MOL3/MOL3/6/MOL6, mature oligodendrocyte subtypes; GWAS, genome-wide association study; SNP, single nucleotide polymorphism; TSS, transcription start site; MAGMA, Multi-marker Analysis of GenoMic Annotation; scDRS, single-cell disease-risk score; MWU, Mann–Whitney U test.

Figure S7 | The oligodendrocyte-lineage risk axis is size-robust, channel-specific, and the S1PR5 detection call survives being resolved one donor at a time.

(A) Rank ladder of the 74 GWAS gene sets on the mature-oligodendrocyte-minus-precursor axis (Δ1) in the four human post-mortem cohorts (GSE277430, GSE118257, GSE301908, GSE231585), shown twice: the rank as computed, and the rank after the regulatory-domain-size trend is removed (colour = rank of 74; the shared rank bar sits above E). Rows are ordered by regulatory-domain size (n_peak descending), which places MS on row 37 and PLATELET on row 56; this is a size ordering, not a rank ordering. The MS rank under each column is 2 / 4 / 7 / 2 as computed and 1 / 1 / 3 / 1 with the size trend removed.

(B) Size-matched tight pool: for each cohort, MS plus the 15 gene sets whose peak-domain size is nearest to MS's, ranked within that pool of 16. MS ranks first in all four cohorts (1/16 in each). The ±15% window used for the sensitivity analysis (40 gene sets) and its ρ(ln n_peak, Δ1) are in Table S7.

(C) Δ1 concordance between the two channels profiled from the same donors (RNA: GSE231585; chromatin: GSE231584 snATAC, ArchR gene scores; the two assays were built separately, so this is a trait-level correspondence, not a paired measurement), 74 gene sets. Spearman ρ falls from 0.722 (raw) to 0.137 (P = 0.24) once both channels are residualised on regulatory-domain size, so the unadjusted cross-channel correlation is dominated by gene-set size. On the plotted axes MS is rank 2 in RNA and rank 36 in the size-matched chromatin series. (In the donor-paired table the same two ranks are 2 → 21; that series is not plotted here and the paired ranks are in Table S7.)

(D) The same comparison in the mouse EAE spinal-cord multiome (GSE250589), where both channels were measured in the same nuclei (14,928 nuclei, 25 libraries, identical cell-type labels, gene sets and seed): ρ = 0.248, P = 0.033; MS ranks 11 in RNA and 36 in chromatin.

(E) The 74-gene-set axis table, 74 gene sets × 21 columns (Table S7) — seven channels (the four human RNA cohorts, human ATAC GSE231584, mouse RNA GSE250589, mouse ATAC GSE250589) × three axes (Δ1, Δ2A, Δ2B); colour = rank of 74 on the shared rank scale; the seven channel strips above the matrix and the three axis strips below it carry the colours keyed in the top band and in the band above E. Rows are ordered by the mean of the four human cohorts' as-computed Δ1 rank (MS on row 2). Note that this is a third row key: Fig 7B orders its rows by the mean size-adjusted rank (MS on row 1), and A orders by domain size (MS on row 37). The MS row is marked in red; AD and SCZ rows in green and blue.

(F) Drug-target × cell-type detection map for GSE118257 (17,799 nuclei; 9 donors from 20 libraries; per-population n printed under each column: OPCs 352, COPs 242, ImOlGs 207, Oligo1 1,129, Oligo2 1,839, Oligo3 775, Oligo4 1,579, Oligo5 1,167, Oligo6 1,484, AST 1,242, MG 796, LYM 423, Endo 836, Peri 697, NEU 5,031 – the four populations that appear in the axis definitions use the same symbols as the delta key and Fig 7C; the OL subtypes and the vascular labels keep the annotation table's own names, because they are a source vocabulary rather than a contrast). Marker diameter = percentage of nuclei in which the target is detected (≥ 1 read; guide circles 1% / 10% / 50%); colour = mean log2(CP10k+1) among detecting nuclei (2.2–4.4 scale, depth-corrected as log2(count × 10,000 / detected genes + 1)); an open circle marks a population with no detecting nucleus. S1PR5 is detected in 6.8–14.9% of differentiated oligodendrocyte states and in COPs (14.9%) but in only 0.85% of OPCs and ≤ 3.6% of any non-lineage population; S1PR1 is endothelial/pericyte-biased (24.8% / 19.9%) and ≤ 1.2% in the lineage; the lymphocyte targets CD52 (0% in all nine lineage columns, 15.8% in immune) and ITGA4 (≤ 1.7% in the lineage, 9.5% in immune) are the specificity controls; NFE2L2 is broadly detected (6.4–38.8%) and is the non-specific reference. MS4A1 is absent from the GSE118257 count matrix, so the fifth row is ITGA4.

(G) The F call resolved per donor. For each of the 27 evaluable donors of the three human cohorts that carry non-lineage populations, the plotted value is that donor's own within-donor contrast: S1PR5 detection percentage in mature oligodendrocytes minus its percentage in the strongest other population of the same donor (percentage points). A population enters a donor's comparison only if that donor contributes ≥ 30 nuclei of it, and a donor is plotted only if its OL pool is ≥ 30 nuclei and at least one comparator survives. The red dashed line is the +1.0 pp margin arm of the rule Fig 7D applies to pooled data, and the tinted band above it is where that arm holds; the grey dotted line is zero, i.e. OL exactly ties the strongest other population. Filled dots meet the full rule (margin ≥ 1.0 pp and OL ≥ 1.5× the comparator, or the comparator is 0): 22 of 27 donors (8/11 GSE301908, 7/8 GSE118257, 7/8 GSE231585). Dot area is the comparator population's nucleus count; the thick bar is the cohort median. All five failures are thin-comparator donors (31–410 nuclei in the population that beat or nearly beat OL), i.e. they are low-information donors rather than counter-examples. GSE277430 is not represented: its released object is an OL-enriched subset with no non-lineage comparator, and its scores were built on normalised counts. Diagnosis is not encoded (Table S7).

Abbreviations: COP = newly formed oligodendrocyte (NFOL); ImOlGs = immature OL/precursor-glial subtypes; Oligo1–6 = the six human mature-OL subtypes of the source study; MG = microglia; Peri = pericyte/smooth muscle; AST = astrocyte; Endo = endothelial; LYM = lymphocytes; NEU = neuronal; OPC = oligodendrocyte precursor cell; EAE = experimental autoimmune encephalomyelitis.

Figure 1

Figure 2

Figure 3

Figure 4

Figure 5

Figure 6

Figure 7

Figure S1

Figure S2

Figure S3, part 1 of 2 (A–D)

Figure S3, part 2 of 2 (E-J)

Figure S4, part 1 of 2

Figure S4, part 2 of 2

Figure S5, part 1 of 3

Figure S5, part 2 of 3

Figure S5, part 3 of 3

Figure S6

Figure S7
