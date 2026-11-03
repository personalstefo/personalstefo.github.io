# Load packages
library(igraph); library(ggraph); library(tidygraph)

load("./data/exercise2.Rdata")

p <- data.frame(interactions)
p_attr <- data.frame (node_attributes)

n <- graph_from_data_frame(p, directed = FALSE, vertices = p_attr)

ggraph(n, layout = "fr") + geom_edge_link(alpha = 0.3) +
  geom_node_point(color = "steelblue") +
  geom_node_text(aes(label = Gene_Id), repel = TRUE, max.overlaps = 1000, size = 2) + theme_void()

ggraph(n, layout = "kk") + geom_edge_link(alpha = 0.3) +
  geom_node_point(color = "steelblue") +
  geom_node_text(aes(label = Gene_Id), repel = TRUE, max.overlaps = 1000, size = 2) + theme_void()

# neighbors() in igraph looks up vertices by the special attribute name. We have Gene_ID, not name. So here we assign all Gene_Ids to the atrribute name instead.
V(n)$name <- V(n)$Gene_Id
tp <- V(n)[name == "TP53"] #V (n) gives us the set of nodes of n, name == "TP53" is a filter on a vertex attribute called name. V(n)[ … ] subsets the vertices using that filter. Result: tp is a vertex sequence (an igraph object, just a selection of vertices that stays tied to a specific graph) containing the single vertex whose name is "TP53".
tp # just to check
tp_neighbors <- neighbors(n, tp, mode = "all") # vertex sequence of neighbors

# What it does:
# neighbors() returns a vertex sequence of nodes directly connected to TP53.
# mode = "all" is fine for undirected graphs (or use "in"/"out" for directed).

length(tp_neighbors) # number of first-order connections
stopifnot(length(tp_neighbors) == degree(n, v = tp, mode = "all")) # This is a line that ChatGPT wrote when I was asking for help. It is useful to think about it a litle. It simply stops R and tells it to throw an error if the length of tp_neighbors that we got, which is the number of first-order connecitons, does not equal the degree of the node tp (degree being the number of incident edges of tp, being the same as the number of neighbors in simple undirected graphs with no self-loops, etc), which is outputted by degree(n, v = tp, mode = "all").

# TRUE only for TP53
V(n)$is_tp53 <- V(n)$name == "TP53"

# "no" for everyone, then "yes" for neighbors of TP53
V(n)$interacts_tp53 <- rep("no", vcount(n))

#V(n)$interacts_tp53 <- …
#V(n) = the vertex sequence (all nodes) of graph n.

#$interacts_tp53 creates/overwrites a vertex attribute called interacts_tp53.

#A vertex attribute is just a vector with one value per node (length must equal vcount(n)), stored on the graph.

#rep("no", vcount(n)) -- vcount(n) returns the number of vertices (nodes) in graph n (419 in this case). rep("no", 419) builds a character vector c("no","no",…,"no") of length 419.

#So we are explicitly creating a value for every node: everyone starts as "no" (i.e., “does not interact with TP53”).

#This is clearer/safer than assigning a single "no" and relying on R’s recycling; it guarantees exact length matching to vcount(n).

#After this line, every node has interacts_tp53 == "no".

#Then you flip the subset of TP53’s neighbors to "yes":

V(n)[tp_neighbors]$interacts_tp53<- "yes"

# quick check
sum(V(n)$is_tp53)              # expect 1
table(V(n)$interacts_tp53)     # "yes" count should match length(tp_neighbors)
all(V(n)[tp_neighbors]$interacts_tp53 == "yes")   # should be TRUE

ggraph(n, layout = "fr") +
  geom_edge_link(alpha = 0.08) +
  geom_node_point(aes(color = interacts_tp53, shape = is_tp53), size = 2.5) +
  scale_color_manual(values = c(no = "grey70", yes = "tomato")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 15)) +
  geom_node_text(aes(label = ifelse(is_tp53 | interacts_tp53 == "yes", name, "")),
                 repel = TRUE, size = 2, max.overlaps = 1000) +
  theme_void()

keep_vs <- c(tp, tp_neighbors) # Concatenates all the neighbors and tp into a single vertex sequence. These are the ones we are keeping.

tp53_ego <- induced_subgraph(n, vids = keep_vs) #the induced_subgraph function takes n, and builds a new graph with only the vertices I decide to keep. These are denoted by the argument vids (stands for vertex identifiers). the object is called ego, because an ego-network (or ego-centric network) is the local subgraph centered on a focal actor (the ego) and their direct connections (alters).


ggraph(tp53_ego, layout = "fr") +
  geom_edge_link(alpha = 0.3) + 
  geom_node_point(aes(color = ifelse(name == "TP53", "TP53", "Partner")), size = 3) +
  scale_color_manual(values = c(TP53 = "tomato", Partner = "steelblue")) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) # label by 'name' (gene symbol) +
  theme_void()

keep_vs <- c(tp, tp_neighbors) # so this is the same across both
drop_vs <- V(n)[!(V(n) %in% keep_vs)] # but now there is this line, which makes it less straightforward compared to induced_subgraph(). this is basically everything in n, but not in keep_vs. 

tp53_ego2 <- delete_vertices(n, drop_vs)
ggraph(tp53_ego2, layout = "fr") +
  geom_edge_link(alpha = 0.3) + 
  geom_node_point(aes(color = ifelse(name == "TP53", "TP53", "Partner")), size = 3) +
  scale_color_manual(values = c(TP53 = "tomato", Partner = "steelblue")) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) # label by 'name' (gene symbol) +
  theme_void()

tp53_ego3 <- make_ego_graph(n, order = 1, nodes = tp, mode = "all")[[1]] # we dont just use the function here because it outputs a list of GRAPHS, so then we use [[1]] to pull the tp-centered ego network out.

ggraph(tp53_ego3, layout = "fr") +
  geom_edge_link(alpha = 0.3) + 
  geom_node_point(aes(color = ifelse(name == "TP53", "TP53", "Partner")), size = 3) +
  scale_color_manual(values = c(TP53 = "tomato", Partner = "steelblue")) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) # label by 'name' (gene symbol) +
  theme_void()

deg <- degree(n, mode = "all") # degree of every vertex
mean(deg) # average degree
transitivity(n, type = "average") # mean of local clustering coefficients over all nodes.(Skip nodes with degree < 2, since clustering coefficient is undefined there.))
diameter(n, directed = FALSE)  # since the graph is undirected

deg <- degree(n)

hist(deg, # the data vector (degrees for all nodes).
     breaks = 50, # Divide the x-axis (degree values) into 50 bins. More bins → more detail; fewer bins → smoother but less detail.
     main   = "Node Degree Distribution", #The title of the plot
     xlab   = "Degree", #Label of the x-axis. Here it represents the degree value (number of neighbors).
     ylab   = "Frequency", # Label of the y-axis. It shows how many nodes fall into each degree bin.
     col    = "lightblue", #Fills the bars of the histogram with light blue color.
     border = "white" #Sets the outline color of the bars to white, making the bars visually cleaner.
     )
tp_deg <- degree(n, v = V(n)[name == "TP53"])
abline(v = tp_deg, col = "red", lwd = 2)  # adds vertical line at TP53's degree. lwd denotes the line width.

clust <- transitivity(n, type = "local")

hist(clust,
     breaks = 50,
     main   = "Local Clustering Coefficient Distribution",
     xlab   = "Clustering Coefficient",
     ylab   = "Frequency",
     col    = "lightgreen", border = "white")
tp_clust <- clust[which(V(n)$name == "TP53")]
abline(v = tp_clust, col = "red", lwd = 2)
tp_clust

# Running Louvain clustering
cl <- cluster_louvain(n, resolution = 2) #We will start with resolution of 1.

# Each vertex gets a cluster membership
membership_vec <- membership(cl)
# ChatGPT explanation:
# membership(): extracts the cluster ID for each vertex from the community object cl.
# Result: membership_vec is an integer vector of length = number of nodes in your graph.
# Example:
# [1] 1 1 2 2 2 3 1 3 ...
# Each entry corresponds to one vertex (in the same order as V(n)), and the number is the cluster label.
# So if membership_vec[10] = 2, that means vertex 10 belongs to cluster 2.
# Cluster labels are arbitrary integers (1, 2, 3, …). They don’t imply ranking or importance, just grouping.

# How many clusters?
length(unique(membership_vec))
table(membership_vec)   # size of each cluster

ggraph(n, layout = "fr") +
  geom_edge_link(alpha = 0.05) +
  geom_node_point(aes(color = factor(membership_vec)), size = 2) +
  geom_node_text(aes(label = ifelse(name == "TP53", name, "")),
                 repel = TRUE, color = "black", size = 3) +
  theme_void()

sub_nodes <- V(n)[membership_vec == 2]
cluster2 <- induced_subgraph(n, vids = sub_nodes)
# Average degree
deg2 <- degree(cluster2)
mean(deg2)

hist(deg2, # the data vector (degrees for all nodes).
     breaks = 50, # Divide the x-axis (degree values) into 50 bins. More bins → more detail; fewer bins → smoother but less detail.
     main   = "Node Degree Distribution of Cluster 2", #The title of the plot
     xlab   = "Degree", #Label of the x-axis. Here it represents the degree value (number of neighbors).
     ylab   = "Frequency", # Label of the y-axis. It shows how many nodes fall into each degree bin.
     col    = "lightblue", #Fills the bars of the histogram with light blue color.
     border = "white" #Sets the outline color of the bars to white, making the bars visually cleaner.
     )

# Average clustering coefficient
transitivity(cluster2, type = "average")

sub_clust <- transitivity(cluster2, type = "local")

hist(sub_clust,
     breaks = 50,
     main   = "Local Clustering Coefficient Distribution of Cluster 2",
     xlab   = "Clustering Coefficient",
     ylab   = "Frequency",
     col    = "lightgreen", border = "white")


# Diameter
diameter(cluster2, directed = FALSE)

ggraph(cluster2, layout = "fr") +
  geom_edge_link(alpha = 0.1) +
  geom_node_point(color = "steelblue", size = 3) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2) +
  theme_void()
