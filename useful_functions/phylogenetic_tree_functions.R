### source this file as follows:
# malik_h_dir <- "/fh/fast/malik_h/"
# if (Sys.info()[["sysname"]]=="Darwin") { malik_h_dir <- "/Volumes/malik_h/" }
# source( paste0(malik_h_dir, "user/jayoung/git_more_repos/Rtest_and_Rnotes/useful_functions/phylogenetic_tree_functions.R") )


###### addInfoToTree - make a treedata object by combining a phylo object with a tibble that has info on the tips
addInfoToTree <- function(tree, info, colnameForTaxonLabels="taxon") {
    ##### get info in same order as taxa in the tree:
    if (! colnameForTaxonLabels %in% colnames(info)) {
        stop("\n\nERROR - there should be a ",colnameForTaxonLabels," column in the info table\n\n")
    }
    ## check all taxa are in the info table
    tipLabelsInInfoTable <- info %>% 
        select(all_of(colnameForTaxonLabels)) %>% 
        deframe()
    missingTaxa <- setdiff(tree$tip.label, tipLabelsInInfoTable)
    if(length(missingTaxa)>0) {
        stop("\n\nERROR - there are taxon labels in the tree that are not in the info table: ",
             paste(missingTaxa, collapse=","),
             "\n\n")
    }
    # now get info
    desiredRows <- match(tree$tip.label, tipLabelsInInfoTable)
    info_treeorder <- info[desiredRows,] 
    # add info to tree
    tree_withInfo <-  left_join(
        tree, 
        info_treeorder ,
        by=c("label"=colnameForTaxonLabels))
    return(tree_withInfo)
}



######### get_mean_dist_to_tips and choose_clades


###### get_mean_dist_to_tips works on a treedata object to get various info about each node . Adds a column to the tbl_tree mean_tip_dist
get_mean_dist_to_tips <- function(tree) {
    my_tbl <- tree %>% 
        as_tibble() 
    my_tbl <- my_tbl %>% 
        mutate(isTip=isTip(tree, .node=my_tbl$node))
    all_tip_ids <- my_tbl %>% 
        filter(isTip) %>% 
        pull(node)
    all_node_dists <- dist.nodes(tree@phylo)
    new_info <- lapply(1:nrow(my_tbl), function(i) {
        if (my_tbl$isTip[i]) {
            output <- tibble(num_desc_tips=1,
                             tot_edgeLen=0,
                             tot_distToTip=0,
                             mean_distToTip=0,
                             all_descendents=NA_character_,
                             tip_descendents=NA_character_)
            return(output)
        } else {
            my_desc <- phytools::getDescendants(tree@phylo, i)
            my_desc_tips <- intersect(my_desc, all_tip_ids)
            
            ## get dists from this node to all descendent tips
            all_dists_to_tips <- all_node_dists[i,my_desc_tips]
            
            subtree <- extract.clade(tree@phylo, my_tbl$node[i])
            output <- tibble(num_desc_tips=length(my_desc_tips),
                             tot_edgeLen=sum(subtree$edge.length),
                             tot_distToTip = sum(all_dists_to_tips),
                             mean_distToTip = mean(all_dists_to_tips),
                             all_descendents=paste(my_desc, collapse=","), 
                             tip_descendents=paste(my_desc_tips, collapse=",")) 
            return(output)
        }
    })  %>% 
        bind_rows() %>% 
        mutate(norm_edgeLen= tot_edgeLen/num_desc_tips)
    my_tbl <- cbind(my_tbl, new_info) %>% 
        as_tibble() %>% 
        relocate(all_descendents, tip_descendents, .after=isTip) 
    return(my_tbl)
}

# 
# ## example data
# x <- "(((Strix_aluco:8.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);"
# tree.owls.lessSymmetric <- read.tree(text= x)
# tree.owls.lessSymmetric <- as.treedata(tree.owls.lessSymmetric)
# 
# my_dists <- get_mean_dist_to_tips(tree.owls.lessSymmetric)
# 
# ## can still extract a tree, but need to do something quite specific to retain branch lengths:
# my_dists %>% as.treedata(branch.length, label)




#### choose_clades function
## uses output of get_mean_dist_to_tips to pick clades with <= a certain total edge length (or mean dist-to-tip, 
choose_clades <- function(dist_tbl, 
                          dist_threshold, 
                          dist_metric = "mean_distToTip", # or could be tot_edgeLen
                          ## what order to check the tips in. Should matter but I may want to check that.
                          randomize_order = TRUE, seed=NULL,
                          assign_internal_nodes=TRUE,
                          quiet=TRUE) {
    
    dist_tbl_tree <- dist_tbl %>% as.treedata(branch.length, label)
    
    ## deal with the case where the distance threshold is higher than anything in the entire tree - we just assing a single clade
    tree_max_dist <- dist_tbl %>% pull(!!dist_metric) %>% max()
    if(dist_threshold >= tree_max_dist) {
        dist_tbl <- dist_tbl %>% 
            mutate(clade="clade_1") %>% 
            relocate(clade, .after=label)
        return(dist_tbl)
    }
    
    all_tip_IDs <- dist_tbl %>% filter(isTip) %>% pull(node)
    if(!quiet) { cat("orig all_tip_IDs ", paste(all_tip_IDs, collapse=","), "\n") }
    ## randomize order
    if(randomize_order) {
        if(!is.null(seed)) {set.seed(seed)}
        all_tip_IDs <- all_tip_IDs[ sample(length(all_tip_IDs), size=length(all_tip_IDs))]
        if(!quiet) { cat("rand all_tip_IDs ", paste(all_tip_IDs, collapse=","), "\n") }
    }
    tip_assignments <- tibble(tip_ID = all_tip_IDs,
                              clade = NA_integer_)
    clade_counter <- 1
    for (each_tip in all_tip_IDs) {
        if(!quiet) { cat("Checking tip ", each_tip, "\n") }
        ## if it's already assigned, we skip it
        assn <- tip_assignments %>% filter(tip_ID==each_tip) %>% pull(clade)
        if(!is.na(assn)) {
            if(!quiet) { cat("    already assigned\n") }
            next
        }
        ## if not, we start looking back at ancestor nodes until we find one whose depth exceeds the treshold
        threshold_exceeded <- FALSE
        this_anc_ID <- each_tip
        while(!threshold_exceeded) {
            ## go to ancestor of current node
            prev_node_checked <- this_anc_ID
            this_anc_ID <- dist_tbl %>% filter(node==this_anc_ID) %>% pull(parent)
            ## get dist_metric (usually mean_distToTip)
            this_dist <- dist_tbl %>% 
                filter(node==this_anc_ID) %>% 
                pull(!! dist_metric )
            ## if we exceed the threshold, assign that clade to all descendent nodes of the PREVIOUS node we checked
            if (this_dist > dist_threshold) {
                if(!quiet) { cat("    node", this_anc_ID, "forms a clade!\n") }
                all_descendents <- phytools::getDescendants(dist_tbl_tree@phylo, prev_node_checked)
                ## only want tips
                all_descendents <- intersect(all_descendents, all_tip_IDs)
                if(!quiet) { cat("       clade members: ", paste(all_descendents, collapse=","), "\n") }
                ## only reassign unassigned descendents: 
                already_assigned_tips <- tip_assignments %>% filter(!is.na(clade)) %>% pull(tip_ID)
                all_descendents <- setdiff(all_descendents, already_assigned_tips)
                if(!quiet) { cat("       clade members without assignment: ", paste(all_descendents, collapse=","), "\n") }
                tip_assignments <- tip_assignments %>% 
                    mutate(clade= case_when(tip_ID %in% all_descendents ~ clade_counter,
                                            TRUE ~ clade))
                clade_counter <- clade_counter + 1
                threshold_exceeded <- TRUE
            }
        }
    }
    ### add that to the input tbl
    dist_tbl <- left_join(dist_tbl, tip_assignments, by=c("node"="tip_ID")) %>% 
        relocate(clade, .after=label) %>% 
        arrange(clade) %>% 
        mutate(clade= case_when( !is.na(clade) ~ paste0("clade_",as.character(clade)),
                                 TRUE ~ NA_character_)) %>%
        mutate(clade=as.factor(clade)) %>% 
        arrange(node)
    
    ### we can also assign clade to any internal node where 100% of the tip descendents have the same clade
    if (assign_internal_nodes) {
        for (each_node in dist_tbl$node) {
            dat <- dist_tbl %>% filter(node==each_node)
            if(dat$isTip) {next}
            if(!is.na(dat$clade)) {next}
            if(!quiet) { cat("## assigning internal node ", each_node, "\n") }
            ## these will be numeric node IDs, including internal nodes
            all_descendents <- phytools::getDescendants(dist_tbl_tree@phylo, each_node) %>% 
                as.integer()
            if(!quiet) { cat("    all_descendents ", paste(all_descendents, collapse = ","), "\n") }
            all_desc_clades <- dist_tbl %>% 
                ## I stumbled on a weird thing that happened because all_descendents is numeric - I think it was assuming I was feeding in a set of colnames, rather than taking me literally.  adding `.env$` is the solution.
                # https://stackoverflow.com/questions/75417188/ambiguity-with-tidyverse-filter-in
                dplyr::filter(node %in% .env$all_descendents) %>%
                filter(!is.na(clade)) %>%
                pull(clade) %>%
                unique()
            if(!quiet) { cat("    all_desc_clades ", paste(all_desc_clades, collapse = ","), "\n") }
            if(length(all_desc_clades)==1) {
                if(!quiet) { cat("        UNIQUE!\n") }
                dist_tbl <- dist_tbl %>% 
                    mutate(clade = case_when(node==each_node ~ all_desc_clades,
                                             TRUE ~ clade))
            }
        }
    }
    
    return(dist_tbl)
}


###### choose_clades_several_thresholds and choose_clades_several_thresholds_report - utility functions to check several thresholds for slicing the tree
choose_clades_several_thresholds <- function(
        dist_tbl, 
        dist_metric = "mean_distToTip", # or could be tot_edgeLen
        distances_to_try = NULL,
        quiet=TRUE) {
    lapply(distances_to_try, function(this_dist) {
        # cat("\n\n####### trying distance ",this_dist, "\n")
        choose_clades(dist_tbl, 
                      dist_threshold=this_dist, 
                      dist_metric=dist_metric,
                      quiet=quiet) %>% 
            mutate(dist_threshold = this_dist)
    }) %>% 
        bind_rows() 
}

choose_clades_several_thresholds_report <- function(
        dist_tbl, 
        dist_metric = "mean_distToTip", # or could be tot_edgeLen
        distances_to_try = NULL,
        quiet=TRUE) {
    
    num_tips <- dist_tbl %>% 
        filter(isTip) %>% 
        pull(node) %>% 
        unique() %>% 
        length()
    
    max_dist <- dist_tbl %>% 
        pull(!!dist_metric) %>% 
        max(na.rm=TRUE)
    
    clade_outputs <- choose_clades_several_thresholds(
        dist_tbl, 
        dist_metric = dist_metric, # or could be tot_edgeLen
        distances_to_try = distances_to_try,
        quiet=quiet)
    
    report <- clade_outputs %>% 
        filter(isTip) %>% 
        group_by(dist_threshold) %>% 
        summarize(num_clades=length(unique(clade))) %>% 
        mutate(tree_num_tips = num_tips) %>% 
        mutate(tree_max_dist = max_dist)
    
    return(report)
}

##### choose_clades demos using example above

# choose_clades_several_thresholds_report(tree_dists,
#                                         distances_to_try=c(10,15,16,17,18,19,20,25,30))

# tree2_tbl_1 <- tree_dists %>%
#     choose_clades(dist_threshold=10)
#
# tree2_tbl_2 <- tree_dists %>%
#     choose_clades(dist_metric = "tot_edgeLen",
#                   dist_threshold=20)
#
# p1 <- tree2_tbl_1 %>%
#     as.treedata(branch.length, label) %>%
#     ggtree(aes(color=clade)) +
#     geom_tiplab(show.legend=FALSE) +
#     geom_treescale(width=10) +
#     labs(title="subtree mean dist to tip < 10")
# p2 <- tree2_tbl_2 %>%
#     as.treedata(branch.length, label) %>%
#     ggtree(aes(color=clade)) +
#     geom_tiplab(show.legend=FALSE) +
#     geom_treescale(width=10) +
#     labs(title="subtree tot length <20")
#
# (p1 + p2) +
#     plot_annotation(title="treeio's sample.nwk tree")