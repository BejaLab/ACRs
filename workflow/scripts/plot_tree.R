library(dplyr)
library(tidyr)
library(treeio)
library(ggtree)
library(readxl)
library(photobiology)
library(ggplot2)
library(ggnewscale)
library(castor)
library(ape)
library(phangorn)

list_file <- "metadata/Channelrhodopsins_Updated_List.xlsx"
tree_file <- "analysis/all/iqtree.treefile.midpoint"

output_file <- "tmp.pdf"

to_treedata <- function(tree) {
    class(tree) <- c("tbl_tree", "tbl_df", "tbl", "data.frame")
    as.treedata(tree)
}

add_hsp <- function(tree, colname) {
    colname <- deparse(substitute(colname))
    treedata <- to_treedata(tree)
    categories <- setNames(tree[[colname]], tree[["label"]]) %>%
        `[`(treedata@phylo$tip.label) %>%
        as.factor
    hsp <- hsp_max_parsimony(treedata@phylo, as.numeric(categories), edge_exponent = 0.1) %>%
        `$`("likelihoods") %>%
        as.data.frame %>%
        setNames(levels(categories)) %>%
        mutate(node = 1:n()) %>%
        gather(value, likelihood, -node) %>%
        filter(likelihood > 0.99) %>%
        setNames(c("node", paste0(colname, "_hsp"), paste0(colname, "_hsp_lh")))
    left_join(tree, hsp, by = "node")
}
get_mrca <- function(phylo, tips) {
    getMRCA(phylo, tips) %>%
        replace(is.null(.), NA)
}
add_mrca <- function(tree, colname) {
    colname <- deparse(substitute(colname))
    treedata <- to_treedata(tree)
    mrca <- mutate(tree, my_column = !!as.name(colname)) %>%
        group_by(my_column) %>%
        mutate(is.tip = label %in% treedata@phylo$tip.label) %>%
        mutate(no_data = all(is.na(my_column))) %>%
        mutate(mrca = get_mrca(treedata@phylo, node[is.tip])) %>%
        mutate(mrca = ifelse(no_data | is.na(mrca), node, mrca)) %>%
        group_by(mrca) %>%
        mutate(enough_tips = sum(is.tip) > 2) %>%
        mutate(ifelse(node == mrca & enough_tips, first(na.omit(my_column)), NA)) %>%
        pull
    tree[[paste0(colname, "_mrca")]] <- mrca
    return(tree)
}

metadata <- read_xlsx(list_file, .name_repair = "universal") %>%
    mutate(Sequence.name = sub(",.+", "", Sequence.name)) %>%
    mutate(Sequence.name = gsub("@", "_", Sequence.name)) %>%
    mutate(Maximum..nm = ifelse(is.na(Action.maximum..nm), Absorption.maximum..nm, Action.maximum..nm)) 

tree <- read.iqtree(tree_file) %>%
    as_tibble %>%
    left_join(metadata, by = c(label = "Sequence.name")) %>%
    mutate(Color = unname(w_length2rgb(Maximum..nm))) %>%
    mutate(Category = case_when(Currents %in% c("no photocurrents", "channel") ~ NA_character_, T ~ gsub("[][]", "", Currents))) %>%
    mutate(Symbol = gsub(",.+", "", Symbol)) %>%
    mutate(Symbol_show = ifelse(Currents.reference == "[Oppermann23](in_prep)", Symbol, NA)) %>%
    add_mrca(ChR.group) %>%
    add_hsp(Category)

cat_colors <- list(
    "anion channel" = "indianred",
    "cation channel" = "deepskyblue",
    "potassium channel" = "purple",
    "channel" = "yellow4"
)
p <- ggtree(to_treedata(tree), aes(color = Category_hsp), layout = "ape") +
    scale_color_manual(values = cat_colors) + new_scale_color() +
    geom_nodepoint(aes(x = branch.x, y = branch.y, subset = !is.na(UFboot) & UFboot >= 95), size = 0.2, color = "#4d4d4dff") +
    geom_nodepoint(aes(x = branch.x, y = branch.y, subset = !is.na(UFboot) & UFboot >= 90 & UFboot < 95), size = 0.2, color = "#b3b3b3ff") +
    geom_tippoint(aes(subset = !is.na(Category) & is.na(Color)), color = "darkgray") +
    geom_tippoint(aes(subset = !is.na(Color), color = Color)) + scale_colour_identity() + new_scale_color() +
    geom_tiplab2(aes(label = Symbol_show), hjust = -0.2) +
    # geom_tiplab2(aes(label = paste(ChR.group, label)), hjust = -0.2, size = 1) +
    geom_treescale(width = 0.5) +
    geom_cladelab(mapping = aes(subset = !is.na(ChR.group_mrca), node = node, label = ChR.group_mrca), offset = -0.1)
    # geom_cladelab(node = 1207, label = "test label", angle = 0, fontsize = 8, hjust = 1)
ggsave(output_file, p)
