##### code for https://github.com/YuLab-SMU/ggtree/issues/705


library(ape)
library(ggtree)
library(tidytree)

# Get example tree data, with tipdata
tr <- rtree(10)

## Make table of information about tips
tree_info <- data.frame(label = paste0("t",1:10),
                        group_character = c(rep("group1", 3),
                                            rep("group2", 4),
                                            rep("group3", 3)))
## add factor column:
tree_info[,"group_factor"] <- factor(tree_info[,"group_character"])
## add integer column:
tree_info[,"group_integer"] <- as.integer(gsub("group","",tree_info[,"group_character"]))
## add TRUE/FALSE column:
tree_info[,"in_group2"] <- tree_info[,"group_character"] == "group2"
tree_info

## Join info and tree together:
tr_plus_info <- left_join(tr, tree_info, by="label")


#### Now plot the tree

# Plot, using `in_group2` column (TRUE/FALSE) to subset geom_tippoint - this works:
ggtree(tr_plus_info) +
    geom_tiplab(aes(color=group_character), offset=0.05) +
    geom_tippoint(aes(subset=in_group2, 
                      color=group_character),
                  size=3, 
                  show.legend = FALSE)


# Plot, using a logical test on the `group_integer` column (group_integer==2) to subset geom_tippoint - this works:
ggtree(tr_plus_info) +
    geom_tiplab(aes(color=group_character), offset=0.05) +
    geom_tippoint(aes(subset= (group_integer==2), 
                      color=group_character),
                  size=3, 
                  show.legend = FALSE)



# Plot, using a logical test on the `group_factor` column (group_factor=="group2") to subset geom_tippoint - this fails:
ggtree(tr_plus_info) +
    geom_tiplab(aes(color=group_character), offset=0.05) +
    geom_tippoint(aes(subset= (group_factor=="group2"), 
                      color=group_character),
                  size=3, 
                  show.legend = FALSE)


# Same problem if we try to subset geom_tiplab instead of geom_tippoint
ggtree(tr_plus_info) +
    geom_tiplab(aes(subset= (group_factor=="group2"), 
                    color=group_character), offset=0.05) +
    geom_tippoint(aes(color=group_character),
                  size=3, 
                  show.legend = FALSE)


# Same problem using the group_character column
ggtree(tr_plus_info) +
    geom_tiplab(aes(color=group_character), offset=0.05) +
    geom_tippoint(aes(subset= (group_character=="group2"), 
                      color=group_character),
                  size=3, 
                  show.legend = FALSE)


# But if I store the desired value in a variable it succeeds
group_of_interest <- "group2"
ggtree(tr_plus_info) +
    geom_tiplab(aes(color=group_character), offset=0.05) +
    geom_tippoint(aes(subset= (group_character==group_of_interest), 
                      color=group_character),
                  size=3, 
                  show.legend = FALSE)



### sessionInfo
sessionInfo()


