# Libraries needed
library(tidyverse)
library(DESeq2)
library(fgsea)
library(msigdbr)

# Data
load(file='./data/exercise5_part1.Rdata')

# 1. Set the sample IDs as rownames of the metadata
rownames(lihcMeta) <- lihcMeta$sample

# 2. Check alignment between metadata and count matrix
all(rownames(lihcMeta) %in% colnames(lichGeneCount))   # should be TRUE
all(rownames(lihcMeta) == colnames(lichGeneCount))    # TRUE if already in same order

# 3. If not same order, fix it:
lichGeneCount <- lichGeneCount[, rownames(lihcMeta)]

# 4. Build the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = lichGeneCount,
  colData   = lihcMeta,
  design    = ~ condition
)

dds

vsd <- rlog(dds)   # variance-stabilizing transform
plotPCA(vsd, intgroup = "condition")

dds <- DESeq(dds)
res <- results(dds)

summary(res)

sum(res$padj < 0.05, na.rm=TRUE)   # number of significant genes

load(file="./data/exercise5_part2.Rdata")

dds_gender <- DESeqDataSetFromMatrix(
  countData = lichGeneCount,
  colData   = lihcMetaFull,
  design    = ~ gender + condition
)

dds_gender <- DESeq(dds_gender)
res_gender <- results(dds_gender, contrast = c("condition", "Primary_Tumor", "Solid_Tissue_Normal"))

summary(res_gender)
sum(res_gender$padj < 0.05, na.rm=TRUE)   # how many genes are significant after adjusting for gender

# 1. Gene sets (Biological Process, Ensembl IDs)
BP_df <- msigdbr(species = "human", category = "C5", subcategory = "BP")
BP_list <- split(BP_df$ensembl_gene, BP_df$gs_name)

# 2. Make ranked stats vector (use Ensembl IDs as names)
myStats <- res_gender[,"stat"]
names(myStats) <- rownames(res_gender)

# 3. Remove version suffix
clean_names <- sub("\\.\\d+$", "", names(myStats))
names(myStats) <- clean_names

# 4. Clean up
myStats <- myStats[!is.na(myStats)]
myStats <- sort(myStats, decreasing = TRUE)

# 5. Run fgsea
fgseaRes <- fgseaMultilevel(pathways = BP_list, stats = myStats,
                            minSize = 15, maxSize = 500, scoreType = "std")

# 6. Results
sum(fgseaRes$padj < 0.05)
head(fgseaRes[order(fgseaRes$padj), ], 10)

# sort significant results
fg_sig <- fgseaRes %>%
  arrange(padj) %>%
  filter(padj < 0.05)

# pick the top pathway
chosen_pathway <- fg_sig$pathway[1]

# plot enrichment
plotEnrichment(BP_list[[chosen_pathway]], myStats) +
  ggplot2::labs(title = paste0("Enrichment: ", chosen_pathway))

# extract first 3 leading-edge genes
leading <- fg_sig$leadingEdge[[1]]
first_genes <- head(leading, 3)
first_genes

stats_abs <- abs(myStats)

fg_nondir <- fgseaMultilevel(
  pathways = BP_list,
  stats = stats_abs,
  minSize = 15,
  maxSize = 500,
  scoreType = "pos"
)

fg_nondir <- fg_nondir[order(fg_nondir$padj), ]
head(fg_nondir, 10)

fg_volcano <- fgseaRes %>%
  mutate(sig = padj < 0.05,
         mlog10padj = -log10(padj))

to_label <- fg_volcano %>% arrange(padj) %>% slice_head(n = 10)

library(ggrepel)
ggplot(fg_volcano, aes(x = NES, y = mlog10padj, color = sig)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_text_repel(data = to_label,
                  aes(label = stringr::str_trunc(pathway, 50)),
                  max.overlaps = 30, size = 3) +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "#1f77b4"),
                     labels = c("n.s.", "< 0.05")) +
  labs(x = "NES",
       y = expression(-log[10](padj)),
       color = "Significant",
       title = "Volcano plot of gene-set enrichment") +
  theme_minimal(base_size = 11)
