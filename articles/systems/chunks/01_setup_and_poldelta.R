library(igraph)
library(ggraph)
library(tidygraph)

hemoglobin <- data.frame(from = c("ALPHA_GLOBIN", "ALPHA_GLOBIN", "BETA_GLOBIN"), to = c("ALPHA_GLOBIN", "BETA_GLOBIN", "BETA_GLOBIN"))
g <- graph_from_data_frame(hemoglobin, directed = FALSE, vertices = NULL)
plot(g)

polymerase_subunit <- data.frame(from = c("DPOD1_HUMAN", "DPOD1_HUMAN", "DPOD1_HUMAN", "DPOD2_HUMAN", "DPOD2_HUMAN", "DPOD3_HUMAN"), to = c("DPOD2_HUMAN", "DPOD3_HUMAN", "DPOD4_HUMAN", "DPOD3_HUMAN", "DPOD4_HUMAN", "DPOD4_HUMAN"))
ps_attr <- data.frame(name = c("DPOD1_HUMAN", "DPOD2_HUMAN", "DPOD3_HUMAN", "DPOD4_HUMAN"), GeneID = c("PolD1", "PolD2", "PolD3", "PolD4"), Catalytic = c("yes", "no", "no", "no"), Description = c(
    "DNA polymerase delta catalytic subunit",
    "DNA polymerase delta subunit 2",
    "DNA polymerase delta subunit 3",
    "DNA polymerase delta subunit 4"
), AA = c(1009, 469, 466, 107))
p <- graph_from_data_frame(polymerase_subunit, directed = FALSE, vertices = ps_attr)
plot(p)
V(p)$Description

ggraph(p) +
    geom_edge_link() +
    geom_node_point()

ggraph(p) + geom_edge_link() + geom_node_point(aes(colour = Catalytic, size = AA)) + geom_node_text(aes(label = name), repel = TRUE)

poldelta_extented_interactions <- data.frame(from = c("DPOD1_HUMAN", "DPOD1_HUMAN", "DPOD1_HUMAN", "DPOD1_HUMAN", "DPOD1_HUMAN", "DPOD1_HUMAN", "DPOD1_HUMAN", "DPOD1_HUMAN", "DPOD1_HUMAN", "DPOD2_HUMAN", "DPOD2_HUMAN", "DPOD2_HUMAN", "DPOD2_HUMAN", "DPOD2_HUMAN", "DPOD2_HUMAN", "DPOD2_HUMAN", "DPOD2_HUMAN", "DPOD3_HUMAN", "DPOD3_HUMAN", "DPOD3_HUMAN", "DPOD3_HUMAN", "DPOD4_HUMAN", "DPOD4_HUMAN", "DPOD4_HUMAN", "DPOD4_HUMAN", "PRI1_HUMAN", "WRIP1_HUMAN"), to = c("DPOD2_HUMAN", "DPOD3_HUMAN", "DPOD4_HUMAN", "S7A6O_HUMAN", "TREX2_HUMAN", "DNA2L_HUMAN", "PRI1_HUMAN", "PRI2_HUMAN", "WRIP1_HUMAN", "DPOD3_HUMAN", "DPOD4_HUMAN", "PDIP2_HUMAN", "BACD1_HUMAN", "WRIP1_HUMAN", "DNA2L_HUMAN", "PRI1_HUMAN", "PRI2_HUMAN", "DPOD4_HUMAN", "DNA2L_HUMAN", "PRI1_HUMAN", "PRI2_HUMAN", "DNA2L_HUMAN", "PRI1_HUMAN", "PRI2_HUMAN", "WRIP1_HUMAN", "PRI2_HUMAN", "WRIP1_HUMAN"), confidence = c(1.00, 1.00, 1.00, 0.18, 1.00, 1.00, 1.00, 1.00, 0.52, 1.00, 1.00, 1.00, 1.00, 0.54, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.57, 1.00, 0.65))
e <- graph_from_data_frame(poldelta_extented_interactions, directed = FALSE, vertices = NULL)


ggraph(e, layout = "kk") + geom_edge_link() + geom_node_point() + geom_node_label(aes(label = name), repel = TRUE)

extended_attr <- data.frame(
    name = c(
        "DPOD1_HUMAN", "DPOD2_HUMAN", "DPOD3_HUMAN", "DPOD4_HUMAN",
        "BACD1_HUMAN", "DNA2L_HUMAN", "PDIP2_HUMAN", "PRI1_HUMAN",
        "PRI2_HUMAN", "S7A6O_HUMAN", "TREX2_HUMAN", "WRIP1_HUMAN"
    ),
    GeneID = c(
        "POLD1", "POLD2", "POLD3", "POLD4",
        "PDIP1", "DNA2", "POLDIP2", "PRIM1",
        "PRIM2", "SLC7A6OS", "TREX2", "WRNIP1"
    ),
    Role_in_replication = c(
        "Polymerase", "Polymerase", "Polymerase", "Polymerase",
        "Uncertain", "Helicase", "Uncertain", "Primase",
        "Primase", "Uncertain", "DNA repair", "DNA repair"
    ),
    Cellular_component = c(
        "Nucleus", "Nucleus", "Nucleus", "Nucleus",
        "Nucleus", "Nucleus; Mitochondrion", "Nucleus", "Nucleus (Nucleoplasm)",
        "Nucleus (Nucleoplasm)", "Nucleus", "Nucleus", "Nucleus"
    ),
    DNA_Binding = c(
        "+", "+", "-", "-",
        "-", "+", "+", "+",
        "+", "-", "+", "-"
    ),
    Description = c(
        "DNA polymerase delta catalytic subunit",
        "DNA polymerase delta subunit 2",
        "DNA polymerase delta subunit 3",
        "DNA polymerase delta subunit 4",
        "Polymerase delta-interacting protein 1",
        "DNA replication helicase/nuclease 2",
        "Polymerase delta-interacting protein 2",
        "DNA primase small subunit",
        "DNA primase large subunit",
        "Probable RNA Pol II localization factor",
        "3′→5′ exonuclease, DNA repair",
        "Werner helicase-interacting protein 1"
    )
)

e <- graph_from_data_frame(poldelta_extented_interactions, directed = FALSE, vertices = extended_attr)

ggraph(e, layout = "kk") +

    geom_edge_link() +

    # NODES
    geom_node_point(
        aes(
            colour = Role_in_replication,
            shape  = DNA_Binding
        ),
        size = 6
    ) +
    geom_node_label(aes(label = GeneID), repel = TRUE, box.padding = 0.4) +

    # Node color scheme
    scale_colour_manual(values = c(
        "Polymerase" = "steelblue",
        "Primase"    = "forestgreen",
        "Helicase"   = "orange",
        "DNA repair" = "purple",
        "Uncertain"  = "gray40"
    ), na.translate = FALSE) +

    # Node shape scheme (+ binds DNA, – does not)
    scale_shape_manual(values = c("+" = 16, "-" = 15), na.translate = FALSE) +

    # Optional: hide edge legends
    guides(edge_colour = "none", edge_alpha = "none", edge_width = "none") +

    theme_void()

ggraph(e, layout = "kk") +
    geom_edge_link(
        aes(
            edge_colour = confidence,
            edge_alpha = confidence,
            edge_width = confidence
        )
    ) +
    scale_edge_colour_gradientn(
        colours = c("plum1", "mediumpurple3", "purple4"),
        values  = c(0, 0.3, 0.9, 1),
        name    = "Confidence" # name for the edge colorbar
    ) +
    guides(edge_colour = guide_edge_colourbar(
        title  = "Confidence",
        ticks  = TRUE,
        breaks = c(0.15, 0.6, 0.95),
        labels = c("0.0–0.3 poor", "0.3–0.9 good", "0.9–1.0 excellent")
    )) +
    scale_edge_alpha(range = c(0.2, 1)) +
    scale_edge_width(range = c(0.3, 1.5)) +

    # NODES
    geom_node_point(aes(colour = Role_in_replication, shape = DNA_Binding), size = 6) +
    geom_node_label(aes(label = GeneID), repel = TRUE, box.padding = 0.4) +
    scale_colour_manual(values = c(
        "Polymerase" = "steelblue",
        "Primase"    = "forestgreen",
        "Helicase"   = "orange",
        "DNA repair" = "purple",
        "Uncertain"  = "gray40"
    ), na.translate = FALSE) +
    scale_shape_manual(values = c("+" = 16, "-" = 15), na.translate = FALSE) +
    theme_void()
