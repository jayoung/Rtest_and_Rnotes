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

