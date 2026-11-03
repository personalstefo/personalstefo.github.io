# Libraries needed
library(igraph)
library(msigdbr)
library(fgsea)
library(ggraph)

# The PPI database we will use is InWeb, created at DTU. Note that it has a confidence score (cs-score) telling us how certain we about each interaction.:
load(file='data/exercise4.Rdata')

seeds <- c("ALDH1A2", "BMP2", "CXADR", "GATA4", "HAS2", "NF1", "NKX2-5", "PITX2", "PKD2", "RXRA", "TBX1", "TBX2", "ZFPM1", "ZFPM2")

# Filter edges by confidence score
inWeb_filtered <- subset(inWeb, inweb_cs >= 0.2)

# Build igraph graph
g <- graph_from_data_frame(d = inWeb_filtered, directed = FALSE)

# Quick check
g

# Step 1: get neighbors (1st-order) of all seeds
neighbors_list <- ego(g, order = 1, nodes = seeds, mode = "all")

# Flatten into a single vector of nodes
neighbor_nodes <- unique(unlist(neighbors_list))

# Extract subgraph (the virtual pulldown)
sub_g <- induced_subgraph(g, vids = neighbor_nodes)

sub_g

# Step 2: Filter the virtual pulldown to remove "sticky" hub proteins
# Function takes a subgraph, the full parent graph, and a cutoff threshold.
filter_virtual_pulldown <- function(subGraph, parentGraph, cutoff) {
  
  # --- Check that both graphs have vertex names (proteins identified by name) ---
  if (is.null(V(subGraph)$name) || is.null(V(parentGraph)$name))
    stop("Both graphs must have vertex 'name' attributes.")
  
  # --- Mode: "all" ensures degrees are calculated without direction ---
  mode <- "all"
  
  # --- Extract node names from the subgraph ---
  vn <- V(subGraph)$name
  
  # --- Degree in the full graph: total number of connections each node has globally ---
  deg_full <- degree(parentGraph, v = vn, mode = mode)
  
  # --- Degree in the subgraph: number of connections each node has inside the pulldown ---
  deg_internal <- degree(subGraph, mode = mode)
  
  # --- Combine results into a data frame for each node ---
  res <- data.frame(
    node          = vn,
    deg_internal  = deg_internal,
    deg_full      = deg_full,
    # Fraction of a node’s total connections that lie within this pulldown
    frac_internal = ifelse(deg_full > 0, deg_internal / deg_full, NA_real_),
    stringsAsFactors = FALSE
  )
  
  # --- Keep only nodes that meet the cutoff for internal fraction ---
  nodes_to_keep <- res$node[!is.na(res$frac_internal) & res$frac_internal >= cutoff]
  
  # --- Build a new subgraph using only the kept nodes ---
  filteredSubGraph <- induced_subgraph(subGraph, vids = V(subGraph)[name %in% nodes_to_keep])
  
  # --- Add the fraction value as a vertex attribute (for later inspection) ---
  V(filteredSubGraph)$frac_internal <- res$frac_internal[match(V(filteredSubGraph)$name, res$node)]
  
  # --- Return the filtered network ---
  filteredSubGraph
}

# Apply filter with different cutoffs
sub_0   <- filter_virtual_pulldown(sub_g, g, cutoff = 0)
sub_02  <- filter_virtual_pulldown(sub_g, g, cutoff = 0.2)
sub_05  <- filter_virtual_pulldown(sub_g, g, cutoff = 0.5)
sub_1   <- filter_virtual_pulldown(sub_g, g, cutoff = 1)

# Summarize number of nodes and edges for each cutoff
sizes <- sapply(
  list(no_cutoff = sub_0, cutoff_0.2 = sub_02, cutoff_0.5 = sub_05, cutoff_1 = sub_1),
  function(x) c(nodes = vcount(x), edges = ecount(x))
)
sizes

# too large for Kamada-Kawai
ggraph(sub_02, layout = "fr") +
  geom_edge_link(alpha = 0.3) + 
  geom_node_point(aes(size = degree(sub_02), color = ifelse(name %in% seeds, "Seed", "Neighbour")), size = 3) +
  scale_color_manual(values = c("Seed" = "tomato", "Neighbour" = "steelblue")) +
  geom_node_text(aes(label = name), repel = FALSE, size = 1) # label by 'name' (gene symbol)

set.seed(22140)  # reproducibility


cl <- cluster_louvain(sub_02, resolution = 0.4)
V(sub_02)$cluster <- membership(cl) # This takes the membership vector (a named integer vector where each element is a node, and the value is the community ID it belongs to) and assigns it as a vertex attribute named "cluster" to the sub_02. 
cluster_sizes <- sort(table(V(sub_02)$cluster), decreasing = TRUE)
cluster_sizes[1:3] # see 3, we will only take the top 2-3 out of the 5 found

# Split into separate subgraphs for each cluster
sub_c1 <- induced_subgraph(sub_02, V(sub_02)$cluster == 1)
sub_c3 <- induced_subgraph(sub_02, V(sub_02)$cluster == 3)
sub_c4 <- induced_subgraph(sub_02, V(sub_02)$cluster == 4)

# Function to plot a single cluster
plot_cluster <- function(graph, title) {
  ggraph(graph, layout = "kk") +
    geom_edge_link(alpha = 0.1, colour = "grey40") +
    geom_node_point(size = 3, aes(color = name %in% seeds)) +
    geom_node_text(aes(label = name), 
                   repel = FALSE, size = 2, color = "black") +
    scale_color_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue"),
                       labels = c("Neighbour", "Seed")) +
    ggtitle(title) +
    theme_void()
}

# Plot clusters separately
plot_cluster(sub_c1, "Cluster 1")
plot_cluster(sub_c3, "Cluster 3")
plot_cluster(sub_c4, "Cluster 4 (largest)")

universe <- V(sub_02)$name
BP_df <- msigdbr(species = "human", category = "C5", subcategory = "BP")

BP_list <- split(BP_df$gene_symbol, BP_df$gs_name)

# Now we need to clean the gene sets, just like in Week 3. 
clean_sets <- function(pathways, universe){
  sets <- lapply(pathways, function(g) intersect(unique(g), universe))
  sets[vapply(sets, length, integer(1)) > 0]
}
BP_list <- clean_sets(BP_list, universe = V(sub_02)$name)

# We know the largest clusters, so now we just take the gene names of the genes in each of those clusters. 

cl1_genes <- V(sub_02)$name[V(sub_02)$cluster == 1]
cl3_genes <- V(sub_02)$name[V(sub_02)$cluster == 3]
cl4_genes <- V(sub_02)$name[V(sub_02)$cluster == 4]

# And run ORA with fora()! 
res_cl1 <- fora(BP_list, genes = cl1_genes, universe = universe)
res_cl3 <- fora(BP_list, genes = cl3_genes, universe = universe)
res_cl4 <- fora(BP_list, genes = cl4_genes, universe = universe)

# Then we write a simple function to get the enrichment

add_enrichment <- function(res, universe, target){
  expected   <- (res$size / length(universe)) * length(target)
  enrichment <- res$overlap / expected
  out <- cbind(res, expected, enrichment)
  out[order(res$padj, res$pval), ]
}

res_cl1 <- add_enrichment(res_cl1, universe, cl1_genes)
res_cl3 <- add_enrichment(res_cl3, universe, cl3_genes)
res_cl4 <- add_enrichment(res_cl4, universe, cl4_genes)

head(res_cl1[, c("pathway","overlap","size","expected","enrichment","pval","padj")], 10)
head(res_cl3[, c("pathway","overlap","size","expected","enrichment","pval","padj")], 10)
head(res_cl4[, c("pathway","overlap","size","expected","enrichment","pval","padj")], 10)

# Betweenness in the largest community (cluster 4)
sub_c4 <- induced_subgraph(sub_02, V(sub_02)$cluster == 4)
bc_c4  <- betweenness(sub_c4, directed = FALSE, normalized = TRUE)

# Top by betweenness
top_gene_c4 <- names(which.max(bc_c4))
top_gene_c4
