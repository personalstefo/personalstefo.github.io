library(tidyverse)
library(fgsea)
library(msigdbr)

# Here I just import the exercise data.
cluster <- read_csv2("data/clusterexGO.csv")
bg_terms <- read_csv2("data/background_countsexGO.csv")
bg_total <- 5500

cluster
bg_terms

# nrow is a base R function that tells the number of rows, in my case it is the number of proteins in the cluster -- 16.
n <- nrow(cluster)
n

# Now we will sum up all the columns. To understand what this does, you have to know the structure of the table. we have a protein X, and then the association of that protein X with three columns GO terms, namely DNA replication (GO:0006260), DNA repair (GO:0006281), Cell cycle (GO:0007049), denoted by 1 (yes, associated) and 0 (no, not associated). So when we sum up all the columns we get the number of proteins in the cluster that are associated with a given GO term. 

observed <- colSums(cluster[, c("DNA_replication", "DNA_repair", "Cell_cycle")])
observed

# bg_terms has the total number of proteins listed that are involved in the three GO categories "DNA replication", "DNA repair" and "Cell cycle" across the entire genome.

# %>% is the pipe operator from the tidyverse. It means take the thing on the left and pass it into the function on the right. 

# mutate() adds a new column to the data frame or overwrites an existing one. 

# So here we’re creating a new column called freq_bg. Inside mutate(), freq_bg = bg_count / bg_total means: divide the number of genes annotated to the term (bg_count) by the total number of background genes (bg_total = 5500). So as a result, for each GO term, we now have the genome-wide frequency.

bg_terms <- bg_terms %>%
  mutate(freq_bg = bg_count / bg_total)

# See the new column.
bg_terms

# Now we do the same -- add a new column, but that column is now the number of expected proteins associated with a given GO term in a randow draw of 16 proteins (n).

bg_terms <- bg_terms %>%
  mutate(expected = freq_bg * n)
bg_terms

# Here we add the observed (colSums) column to the table matched by name of GO term. 
bg_terms <- bg_terms %>%
  mutate(observed = observed[term])
bg_terms

# So enrichment is simply the number of observed divided by the number of expected, which we calculated earlier. Reminder: expected is a decimal, but not a proportion!

bg_terms <- bg_terms %>%
  mutate(enrichment = observed/expected)
bg_terms

# Observed associated with DNA replication
a <- observed["DNA_replication"]
# Number of proteins in cluster NOT associated with DNA replication
b <- n - a
# "c" is the number of DNA replication genes OUTSIDE the cluster

# Step 1:bg_terms$bg_count[] looks up how many DNA replication genes exist in the whole genome (in this case: 96 total). On the operator $: bg_terms$bg_count means: “From the data frame bg_terms, give me the column named bg_count.”. when we specify bg_terms$term == "DNA_Replication", we get the bg_count of dna_replication.

# Step 2: subtract "a", the ones we already counted inside the cluster. That leaves only the ones not in the cluster.
c <- bg_terms$bg_count[bg_terms$term == "DNA_replication"] - a

# here the login is the same, but it is the number of genes that in the whole genome are NOT associated with dna replication minus the genes in the cluster that are not associated with DNA replication.

d <- (bg_total - bg_terms$bg_count[bg_terms$term == "DNA_replication"]) - b

# Here we create the table mentioned above.
m <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
            dimnames = list(c("has_term","no_term"),
                            c("in_cluster","not_in_cluster")))

m

# alternative = "greater" sets the direction of the hypothesis test.
# "greater" = test if the GO term is over-represented in the cluster.
# (Other options: "less" tests under-representation; "two.sided" tests either way.)

fisher.test(m, alternative = "greater")

# DNA repair
a <- observed["DNA_repair"]
b <- n - a
c <- bg_terms$bg_count[bg_terms$term == "DNA_repair"] - a
d <- (bg_total - bg_terms$bg_count[bg_terms$term == "DNA_repair"]) - b

m_repair <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
            dimnames = list(c("has_term","no_term"),
                            c("in_cluster","not_in_cluster")))

fisher.test(m_repair, alternative = "greater")

# Cell cycle
a <- observed["Cell_cycle"]
b <- n - a
c <- bg_terms$bg_count[bg_terms$term == "Cell_cycle"] - a
d <- (bg_total - bg_terms$bg_count[bg_terms$term == "Cell_cycle"]) - b

m_cycle <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
            dimnames = list(c("has_term","no_term"),
                            c("in_cluster","not_in_cluster")))

fisher.test(m_cycle, alternative = "greater")

load("data/exercise3_part1.Rdata")

load("data/exercise3_part2.Rdata")

# From now on, I will personally need a lot of help from ChatGPT, but I will try to decipher everything as well.

# This command retrieves biological process gene sets for yeast. The "C5" category is the msigdb annotation for the three gene ontologies. The subcategory is either "BP" (biological process), "MF" (molecular function), or "CC" (cellular component).

# Here we take BP (biological process)

BP_df = msigdbr(species = "S. cerevisiae", category = "C5", subcategory = "BP")

# Inspect structure and a few rows
dim(BP_df)          # how many rows & columns?
names(BP_df)        # column names available
head(BP_df[, c("gs_name","ensembl_gene")], 5)   # key columns we’ll use

BP_list = split(x = BP_df$ensembl_gene, f = BP_df$gs_name)

# Inspect to see what one looks like
length(BP_list)          # number of gene sets
names(BP_list)[1:5]       # first few set names
head(BP_list[[1]], 10)    # genes inside first set

# Universe: all genes that are "possible draws" for this analysis.
# We take them from the node table's ORF IDs so they match msigdbr's ORF IDs.
universe <- unique(node_attributes$ID)   # unique() removes duplicates

# Then clean BP sets to that universe
clean_sets <- function(pathways, universe){
  sets <- lapply(pathways, function(g) intersect(unique(g), universe))
  sets[vapply(sets, length, integer(1)) > 0]
}
BP_list <- clean_sets(BP_list, universe)

# Quick sanity peek:
length(universe)     # How big is the universe? (a count)
head(universe, 5)    # What do these IDs look like? Should be Yxxx... style


# Target: the actual genes we're testing for enrichment (Cluster #1).
# Boolean filter: keep rows where the 'cluster' column equals the label "cluster1",
# then pull their ORF IDs and drop duplicates.
target <- unique(node_attributes$ID[node_attributes$cluster == "cluster1"])

# Quick sanity peek:
length(target)       # How many genes are in Cluster #1?
head(target, 5)      # Do these look like ORF IDs too?

# Sanity checks that prevent nonsense results:
sum(is.na(target))                    # should be 0 (no missing IDs in target)
all(target %in% universe)             # should be TRUE (every target in the universe)

# ⚠️ Problem — one of your target entries is NA. That means at least one row in cluster1 had a missing ID.

# Removing it: 

target <- target[!is.na(target)]

# Rechecking: 

length(target)
all(target %in% universe)

# Run over-representation for BP
res_bp <- fora(pathways = BP_list,   # named list: set -> vector of ORF IDs
               genes    = target,    # your cluster-1 ORF IDs
               universe = universe)  # all possible ORF IDs (background)

# See what columns we got back
names(res_bp)        # should include: pathway, pval, padj, size, overlap
head(res_bp, 5)      # quick peek at the first rows

add_enrichment <- function(res, universe, target){
  expected   <- (res$size / length(universe)) * length(target)  # random-draw expectation
  enrichment <- res$overlap / expected                          # effect size (observed/expected)
  out <- cbind(res, expected = expected, enrichment = enrichment)
  out[order(out$padj, out$pval), ]                              # rank by adjusted p-value
}

res_bp <- add_enrichment(res_bp, universe, target)

# Show top-10 BP terms for your report
top10_bp <- head(res_bp[, c("pathway","overlap","size","expected","enrichment","pval","padj")], 10)
top10_bp

# Get MF & CC from msigdbr
MF_df <- msigdbr(species = "S. cerevisiae", category = "C5", subcategory = "MF")
CC_df <- msigdbr(species = "S. cerevisiae", category = "C5", subcategory = "CC")

# Convert to lists
MF_list <- split(MF_df$ensembl_gene, MF_df$gs_name)
CC_list <- split(CC_df$ensembl_gene, CC_df$gs_name)

# Clean to your universe
MF_list <- clean_sets(MF_list, universe)
CC_list <- clean_sets(CC_list, universe)

# Run fora
res_mf <- fora(MF_list, target, universe)
res_cc <- fora(CC_list, target, universe)

# Add enrichment and take top 10
res_mf  <- add_enrichment(res_mf, universe, target)
res_cc  <- add_enrichment(res_cc, universe, target)

top10_mf <- head(res_mf[, c("pathway","overlap","size","expected","enrichment","pval","padj")], 10)
top10_cc <- head(res_cc[, c("pathway","overlap","size","expected","enrichment","pval","padj")], 10)

top10_mf
top10_cc
