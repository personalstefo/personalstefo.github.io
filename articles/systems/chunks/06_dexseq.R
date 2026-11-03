#Loading environment and data

library(tidyverse)
library(DEXSeq)
library(fgsea)
library(msigdbr)
load(file='./data/exercise6.Rdata')

# Number of isoforms and genes we are analyzing
n_isoforms_input <- length(unique(txInfoSubset$transcript_id))
n_genes_input <- length(unique(txInfoSubset$gene_id))
n_isoforms_input
n_genes_input

dxdSubset <- DEXSeqDataSet(
countData = lichIsoCountSubset,
sampleData = lihcMetaFull,
design= ~ sample + exon + condition:exon + gender:exon,
featureID = txInfoSubset$transcript_id,
groupID = txInfoSubset$gene_id
)

dxdSubset <- estimateSizeFactors(dxdSubset)
dxdSubset <- estimateDispersions(dxdSubset)

dxdSubset <- testForDEU(
dxdSubset,
fullModel = ~ sample + exon + condition:exon + gender:exon,
reducedModel = ~ sample + exon + gender:exon
)

dxdSubset <- estimateExonFoldChanges(dxdSubset, fitExpToVar = "condition")
dxdResult <- DEXSeqResults(dxdSubset)
alpha <- 0.05

# Significant isoforms (rows)
sig_iso <- dxdResult[!is.na(dxdResult$padj) & dxdResult$padj < alpha, ]
n_sig_isoforms <- nrow(sig_iso)
n_sig_isoforms

n_sig_genes <- length(unique(sig_iso$groupID))
n_sig_genes

dxdResult = DEXSeq(dxdSubset)
class(dxdResult)

# Load full analysis
load(file='./data/exercise6_intermediary_results.Rdata')

dexseq_df <- as.data.frame(dexseqResult)
dexseq_df <- dexseq_df[order(dexseq_df$padj), ]
top5_genes <- dexseq_df$groupID[!is.na(dexseq_df$padj)][1:5]
top5_genes

# GOAL:
# Collapse isoform-level DEXSeq results to one row per gene by:
#   (a) picking the most significant isoform per gene (smallest RAW p-value),
#   (b) then recomputing FDR across genes (new hypothesis family).
# WHY:
#   The original DEXSeq 'padj' adjusts across ALL isoforms/features. After collapsing to
#   one p-value per gene, that 'padj' is no longer valid; we must re-adjust at the gene level.

dexseqResultGene <-
  dexseqResult %>%
  as.data.frame() %>%        # Convert S4 results object to a plain data.frame
  as_tibble() %>%            # and then to a tibble for friendly dplyr usage

  # 1) Keep only tested features (rows where the test statistic exists).
  #    NAs often arise from zero counts, filtering, or non-convergence.
  filter(!is.na(stat)) %>%

  # 2) Select and rename columns we need downstream.
  #    - groupID   : gene identifier as supplied when building the DEXSeqDataSet (here: Ensembl gene id)
  #    - featureID : feature identifier (here: transcript/isoform id, because we set featureID = transcript_id)
  #    - stat      : test statistic for the DEU test on this feature
  #    - log2fold_Primary_Tumor_Solid_Tissue_Normal : feature-level log2FC for the 'condition' contrast
  #    - pvalue    : RAW p-value for this feature (used to pick per-gene "best" isoform)
  #    - padj      : feature-level FDR across ALL features (not valid after collapsing to genes)
  dplyr::select(
    gene_id    = groupID,
    isoform_id = featureID,
    stat,
    log2FC     = log2fold_Primary_Tumor_Solid_Tissue_Normal,
    pvalue,
    padj
  ) %>%

  # 3) Group by gene to operate within each gene's set of isoforms/features.
  group_by(gene_id) %>%

  # 4) Within each gene, keep the single isoform with the smallest RAW p-value.
  #    WHY raw p (not feature-level padj)? Because we are about to define a NEW testing family
  #    (one p per gene) and must feed those raw p-values into a fresh multiple-testing correction.
  slice_min(
    order_by = pvalue,
    n = 1,
    with_ties = FALSE  # if two isoforms tie exactly, pick one arbitrarily (see tie-breaking alternative below)
  ) %>%

  # 5) Ungroup so the next mutate acts on the entire tibble.
  ungroup() %>%

  # 6) Recompute FDR across genes (new hypothesis family = "one p per gene").
  #    This is statistically correct and typically more powerful (fewer tests than "all isoforms").
  #    Note: method = "fdr" is an alias for Benjamini–Hochberg ("BH").
  mutate(padj = p.adjust(pvalue, method = "fdr"))


# --------------------------- OPTIONAL: nicer tie-breaking ---------------------------
# For deterministic tie-breaking instead of with_ties = FALSE, we can
# order by p-value first and then by absolute log2FC (or by coverage) and take the first row.
# This can be helpful if multiple isoforms have identical p-values for a gene.

# dexseqResultGene_tiebroken <-
#   dexseqResult %>%
#   as.data.frame() %>%
#   as_tibble() %>%
#   filter(!is.na(stat)) %>%
#   dplyr::select(
#     gene_id    = groupID,
#     isoform_id = featureID,
#     stat,
#     log2FC     = log2fold_Primary_Tumor_Solid_Tissue_Normal,
#     pvalue,
#     padj
#   ) %>%
#   group_by(gene_id) %>%
#   arrange(pvalue, desc(abs(log2FC)), .by_group = TRUE) %>%  # tie-break by larger |log2FC|
#   slice(1) %>%
#   ungroup() %>%
#   mutate(padj = p.adjust(pvalue, method = "fdr"))


# --------------------------- QUICK SANITY CHECKS (optional) ---------------------------
# # All kept rows should have p-values:
# stopifnot(!any(is.na(dexseqResultGene$pvalue)))
#
# # Exactly one row per gene:
# stopifnot(dplyr::n_distinct(dexseqResultGene$gene_id) == nrow(dexseqResultGene))
#
# # Inspect gene-level FDR distribution:
# summary(dexseqResultGene$padj)


# --------------------------- COMMON GOTCHAS ---------------------------
# • ID versions: If the gene IDs include version suffixes (e.g., "ENSG... .12"), we need to keep them consistent across any joins with DESeq2 outputs (strip or keep versions on both sides).
# • Low coverage: Very lowly expressed features can yield noisy min-p selections; consider prefiltering by coverage/baseMean if you see suspicious genes.
# • Non-converged rows: The earlier filter(!is.na(stat)) drops non-converged/untested features, which is intended before ranking and FDR.

###############################################
# STEP 1 — Choose a cutoff for "significant"
###############################################

# Here we decide that any adjusted p-value (padj) smaller than 0.05
# means that the gene is significantly different.
# 0.05 is a common threshold — it represents a 5% false discovery rate.
alpha <- 0.05

###############################################
# STEP 2 — (Optional) Clean gene IDs
###############################################

# Sometimes gene IDs have version numbers at the end (like ENSG000001234.12).
# DESeq2 and DEXSeq might handle these differently, so we make sure both use the same format.
# This function removes everything after a dot and the digits following it.
strip_ver <- function(x) sub("\\.\\d+$", "", x)

###############################################
# STEP 3 — Prepare the DEXSeq results
###############################################

# The dexseqResultGene table already has one row per gene (after collapsing isoforms).
# We add a TRUE/FALSE column saying whether each gene is significant.

dex_gene <- dexseqResultGene %>%
  mutate(
    gene_id = strip_ver(gene_id),   # make sure gene IDs are version-free and comparable
    dex_padj = padj,                # keep padj for clarity
    dex_sig  = !is.na(padj) & padj < alpha  # TRUE if padj < 0.05 and not missing
  ) %>%
  dplyr::select(gene_id, dex_padj, dex_sig) # keep only relevant columns
dex_gene

###############################################
# STEP 4 — Prepare the DESeq2 results
###############################################

# Do exactly the same thing for the DESeq2 gene-level results.
# Again, add a TRUE/FALSE flag for significance.

deseq_gene <- deseq2Results %>%
  mutate(
    gene_id = strip_ver(gene_id),   # remove version suffix
    deseq_padj = padj,              # keep padj column
    deseq_sig  = !is.na(padj) & padj < alpha  # mark significant genes
  ) %>%
  dplyr::select(gene_id, deseq_padj, deseq_sig)
deseq_gene

###############################################
# STEP 5 — Define the universe of comparable genes
###############################################

# Both DEXSeq and DESeq2 may have tested slightly different sets of genes.
# We only want to compare the ones that appear in BOTH.
# This gives us the "gene universe" for the overlap test.
univ <- intersect(dex_gene$gene_id, deseq_gene$gene_id)

# Sanity check: how many genes are shared between the two analyses?
length(univ)

###############################################
# STEP 6 — Combine both results side-by-side
###############################################

# Merge the two tables by gene_id, so that each gene now has:
# - its DEXSeq significance (dex_sig)
# - its DESeq2 significance (deseq_sig)

both <- dex_gene %>%
  filter(gene_id %in% univ) %>%      # keep only shared genes
  inner_join(deseq_gene, by = "gene_id")  # join DESeq2 results by gene_id

# Quick look to make sure it worked
head(both)

###############################################
# STEP 7 — Count overlaps (the "both TRUE" genes)
###############################################

# This is the number of genes significant in both methods.
n_overlap <- sum(both$dex_sig & both$deseq_sig)

# Also useful to know total counts for context.
n_dex_sig   <- sum(both$dex_sig)
n_deseq_sig <- sum(both$deseq_sig)
n_universe  <- nrow(both)

# Print them
cat("Universe:", n_universe, "\n")
cat("DEXSeq significant:", n_dex_sig, "\n")
cat("DESeq2 significant:", n_deseq_sig, "\n")
cat("Overlap:", n_overlap, "\n")


###############################################
# STEP 8 — Build the 2×2 contingency table
###############################################

# This table shows every possible combination:
# - Neither method found it
# - Only DESeq2 found it
# - Only DEXSeq found it
# - Both found it

contig_tab <- table(DEXSeq = both$dex_sig, DESeq2 = both$deseq_sig)
contig_tab

###############################################
# STEP 9 — Run Fisher’s exact test (ORA)
###############################################

# Fisher’s Exact Test asks:
# "Is the overlap (TRUE + TRUE) larger than we would expect by random chance?"

ft <- fisher.test(contig_tab, alternative = "greater")

# Print the result
ft

###############################################
# STEP 10 — Make a tidy summary table
###############################################

tidy_summary <- tibble(
  universe     = n_universe,
  dex_sig      = n_dex_sig,
  deseq_sig    = n_deseq_sig,
  overlap      = n_overlap,
  odds_ratio   = unname(ft$estimate),
  p_value      = ft$p.value
)
tidy_summary

set.seed(228)
# 1. Gene sets (Biological Process, Ensembl IDs)
BP_df <- msigdbr(species = "human", category = "C5", subcategory = "BP")
BP_list <- split(BP_df$ensembl_gene, BP_df$gs_name)

# 2. Create ranked gene list for DESeq2
# Rank all genes by their test statistic ('stat'), which represents how strongly
# they are differentially expressed (positive = up, negative = down).

stats_deseq <- deseq2Results$stat
names(stats_deseq) <- sub("\\.\\d+$", "", deseq2Results$gene_id)  # clean Ensembl IDs
stats_deseq <- stats_deseq[!is.na(stats_deseq)]  # drop NAs
stats_deseq <- sort(stats_deseq, decreasing = TRUE)  # rank from high to low


# 3. Create ranked gene list for DEXSeq (gene-level)
# Same idea, but now using DEXSeq gene-level results (from the collapsed table).

stats_dex <- dexseqResultGene$stat
names(stats_dex) <- sub("\\.\\d+$", "", dexseqResultGene$gene_id)
stats_dex <- stats_dex[!is.na(stats_dex)]
stats_dex <- sort(stats_dex, decreasing = TRUE)

# 4. Run GSEA for each analysis
# fgseaMultilevel tests each pathway to see if its genes cluster at the top or bottom
# of our ranked list more than expected by chance.

fg_deseq <- fgseaMultilevel(pathways = BP_list, stats = stats_deseq,
                            minSize = 15, maxSize = 500, scoreType = "std")

fg_dex <- fgseaMultilevel(pathways = BP_list, stats = stats_dex,
                          minSize = 15, maxSize = 500, scoreType = "std")

# 5. Count how many significant pathways each method found
alpha <- 0.05
n_sig_de <- sum(fg_deseq$padj < alpha, na.rm = TRUE)
n_sig_dx <- sum(fg_dex$padj < alpha, na.rm = TRUE)

cat("Significant pathways (DESeq2):", n_sig_de, "\n")
cat("Significant pathways (DEXSeq):", n_sig_dx, "\n")

shared <- intersect(
  fg_deseq$pathway[fg_deseq$padj < 0.05],
  fg_dex$pathway[fg_dex$padj < 0.05]
)
cat("Shared pathways:", length(shared), "\n")

############################################################
# GOAL
# Compare pathway-level enrichment between DESeq2 (gene-level)
# and DEXSeq (isoform-usage) by plotting their NES scores
# against each other. Color points by which analysis called
# the pathway significant (DEXSeq only, DESeq2 only, Both, None).
#
# INPUTS WE ALREADY HAVE:
#   fg_deseq : fgseaMultilevel() result for DESeq2 stats
#   fg_dex   : fgseaMultilevel() result for DEXSeq  stats
#   alpha    : significance cutoff (e.g., 0.05)
#
# OUTPUT:
#   1) A merged data.frame `cmp` with NES + padj from both analyses
#   2) A scatterplot NES(DESeq2) vs NES(DEXSeq) colored by significance class
############################################################

library(dplyr)
library(ggplot2)

# 1) Keep just the columns we need from each fgsea result and rename them
#    so they won't collide when we join the two tables later.
#    - pathway : the gene set name (join key)
#    - NES     : normalized enrichment score (effect size; up/down and magnitude)
#    - padj    : FDR-adjusted p-value for the gene set (significance)
#
#    `transmute()` creates a new data.frame with only the columns we list.
#    We also rename NES/padj to method-specific names (NES_de/padj_de, NES_dx/padj_dx).
de_df <- fg_deseq %>%
  transmute(
    pathway,
    NES_de  = NES,
    padj_de = padj
  )

dx_df <- fg_dex %>%
  transmute(
    pathway,
    NES_dx  = NES,
    padj_dx = padj
  )

# 2) Merge (inner join) on 'pathway' so we keep only gene sets present in BOTH
#    results. This is necessary because we want to plot NES_de vs NES_dx for
#    the same pathway on each row.
#
#    After the join, add a "sig_cat" (significance category) column:
#      - "Both"         : padj < alpha in DESeq2 AND in DEXSeq
#      - "DESeq2 only"  : padj < alpha in DESeq2, but not in DEXSeq
#      - "DEXSeq only"  : padj < alpha in DEXSeq, but not in DESeq2
#      - "None"         : padj >= alpha in both
#
#    Finally, we make sig_cat an ordered factor so the legend uses a sensible order.
cmp <- inner_join(de_df, dx_df, by = "pathway") %>%
  # fgsea can leave a few rows with NA NES/padj (e.g., "unbalanced" sets).
  # Drop NAs so plotting/statistics don't choke.
  filter(!is.na(NES_de), !is.na(NES_dx), !is.na(padj_de), !is.na(padj_dx)) %>%
  mutate(
    sig_cat = case_when(
      padj_de < alpha & padj_dx < alpha ~ "Both",
      padj_de < alpha & padj_dx >= alpha ~ "DESeq2 only",
      padj_de >= alpha & padj_dx < alpha ~ "DEXSeq only",
      TRUE ~ "None"
    ),
    sig_cat = factor(sig_cat, levels = c("Both", "DESeq2 only", "DEXSeq only", "None"))
  )

# 3) (Optional but useful) Compute correlations so we can report a number.
#    Spearman = monotonic rank correlation (robust to outliers, scale-free).
#    Pearson  = linear correlation.
spearman_cor <- cor(cmp$NES_de, cmp$NES_dx, method = "spearman", use = "complete.obs")
pearson_cor  <- cor(cmp$NES_de, cmp$NES_dx, method = "pearson",  use = "complete.obs")

# 4) Make the scatterplot:
#    - x = DESeq2 NES; y = DEXSeq NES
#    - color by significance category
#    - add a linear trend line (geom_smooth with method = "lm") WITHOUT a gray CI ribbon
#      to avoid clutter; the line just indicates the direction/strength visually.
p <- ggplot(cmp, aes(x = NES_de, y = NES_dx, color = sig_cat)) +
  geom_point(alpha = 0.8) +                 # points for each pathway
  geom_smooth(method = "lm", se = FALSE) +  # dashed trend line (no CI)
  labs(
    x = "NES (DESeq2 fgsea)",               # clear axis labels
    y = "NES (DEXSeq fgsea)",
    color = "Significance",
    title = "NES comparison: DESeq2 vs DEXSeq (fgseaMultilevel)",
    subtitle = paste0("Spearman r = ", round(spearman_cor, 3),
                      " | Pearson r = ", round(pearson_cor, 3),
                      " | alpha = ", alpha)
  ) +
  theme_minimal()

# 5) Print the plot
p

# 6) (Optional) Peek at how many sets fall in each category;
table(cmp$sig_cat)
