## https://github.com/YuLab-SMU/treeio/issues/120

## I also tried an older version treeio_1.20.2 (on a different computer) and this time, the tip.label and node.labels are retained. So something went wrong between versions 1.20.2 and version 1.26.0.

# It is NOT a problem with ape::root. 
# It is NOT a problem on my home Mac that has an older version of treeio (treeio_1.20.2).  I think I could downgrade the treeio version I'm using. e.g. to install an old package version:
 ## CRAN
    # packageurl <- "http://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_0.9.1.tar.gz"
    # install.packages(packageurl, repos=NULL, type="source")
 ## bioconductor:
    # https://stackoverflow.com/questions/49487171/r-how-to-install-a-specified-version-of-a-bioconductor-package


# install_github()

# Here's how I view a function: `treeio:::root.treedata`

### In treeio_1.26.0 the `treeio:::root.treedata` command gives me an error:
# Error: object 'root.treedata' not found
### and this command indicates that the method is 'hidden' (the asterisk)
# methods ("root")
# [1] root.multiPhylo* root.phylo*      root.treedata*  


### In treeio_1.20.2 the function looks like this and works as expected.  
# treeio:::root.treedata
# function (phy, outgroup, node = NULL, edgelabel = TRUE, ...)
# {
#     if (!missing(outgroup) && is.character(outgroup)) {
#         outgroup <- match(outgroup, phy@phylo$tip.label)
#     }
#     if (!edgelabel) {
#         message("The use of this method may cause some node data to become incorrect (e.g. bootstrap values) if 'edgelabel' is FALSE.")
#     }
#     res <- build_new_labels(tree = phy)
#     tree <- res$tree
#     node2oldnewlab <- res$node2old_new_lab
#     re_tree <- root(tree, outgroup = outgroup, node = node, edgelabel = edgelabel,
#                     ...)
#     node_map <- old_new_node_mapping(tree, re_tree)
#     n.tips <- Ntip(re_tree)
#     phy@phylo <- build_new_tree(tree = re_tree, node2old_new_lab = node2oldnewlab)
#     update_data <- function(data, node_map) {
#         cn <- colnames(data)
#         cn <- cn[cn != "node"]
#         data <- dplyr::inner_join(data, node_map, by = c(node = "old")) |>
#             dplyr::select(c("new", cn)) |> dplyr::rename(node = .data$new)
#         root <- data$node == (n.tips + 1)
#         data[root, ] <- NA
#         data[root, "node"] <- n.tips + 1
#         return(data)
#     }
#     if (nrow(phy@data) > 0) {
#         phy@data <- update_data(phy@data, node_map)
#     }
#     if (nrow(phy@extraInfo) > 0) {
#         phy@extraInfo <- update_data(phy@extraInfo, node_map)
#     }
#     return(phy)
# }
# <bytecode: 0x7ff2a5647388>
#     <environment: namespace:treeio>
#     


# 
# xxx maybe I don't want to sync to git until I have homing_endos_5_analyzeBaliPhyOutput.Rmd running like it was before?



################### CODE


## load package, get data, 
library(treeio)
data(bird.orders, package="ape")
# add fake node labels
bird.orders$node.label <- paste("node",1:Nnode(bird.orders),sep="")
# convert to treedata
bird.orders_treedata <- as.treedata(bird.orders)

# confirm that there are meaningful tip and node labels
bird.orders_treedata@phylo$tip.label
#  [1] "Struthioniformes" "Tinamiformes"     "Craciformes"      "Galliformes"      "Anseriformes"    
#  [6] "Turniciformes"    "Piciformes"       "Galbuliformes"    "Bucerotiformes"   "Upupiformes"     
# [11] "Trogoniformes"    "Coraciiformes"    "Coliiformes"      "Cuculiformes"     "Psittaciformes"  
# [16] "Apodiformes"      "Trochiliformes"   "Musophagiformes"  "Strigiformes"     "Columbiformes"   
# [21] "Gruiformes"       "Ciconiiformes"    "Passeriformes"   
bird.orders_treedata@phylo$node.label
#  [1] "node1"  "node2"  "node3"  "node4"  "node5"  "node6"  "node7"  "node8"  "node9"  "node10" "node11"
# [12] "node12" "node13" "node14" "node15" "node16" "node17" "node18" "node19" "node20" "node21" "node22"

##### reroot
bird.orders_treedata_rerooted <- root(bird.orders_treedata, "Galliformes")

# we no longer have those meaningful tip and node labels
bird.orders_treedata_rerooted@phylo$tip.label
#  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20"
# [21] "21" "22" "23"

bird.orders_treedata_rerooted@phylo$node.label
#  [1] "24" "25" "26" "27" "28" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "39" "40" "41" "42" "43"
# [21] "44"


sessionInfo()
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.6.1
# 
# Matrix products: default
# BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# locale:
#     [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Los_Angeles
# tzcode source: internal
# 
# attached base packages:
#     [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#     [1] treeio_1.26.0
# 
# loaded via a namespace (and not attached):
#     [1] vctrs_0.6.4       nlme_3.1-163      cli_3.6.1         rlang_1.1.2       purrr_1.0.2      
# [6] generics_0.1.3    jsonlite_1.8.7    glue_1.6.2        fansi_1.0.5       grid_4.3.2       
# [11] tibble_3.2.1      fastmap_1.1.1     ape_5.7-1         lifecycle_1.0.4   memoise_2.0.1    
# [16] compiler_4.3.2    dplyr_1.1.4       fs_1.6.3          Rcpp_1.0.11       pkgconfig_2.0.3  
# [21] tidytree_0.4.5    tidyr_1.3.0       lattice_0.21-9    digest_0.6.33     R6_2.5.1         
# [26] tidyselect_1.2.0  utf8_1.2.4        pillar_1.9.0      parallel_4.3.2    magrittr_2.0.3   
# [31] withr_2.5.2       tools_4.3.2       lazyeval_0.2.2    cachem_1.0.8      yulab.utils_0.1.0


####### test other versions:


#### treeio_1.24.3 has the problem (bioc 3.17 default version)
# install.packages("https://bioconductor.org/packages/3.17/bioc/src/contrib/treeio_1.24.3.tar.gz", repos=NULL, type="source")

#### treeio_1.22.0 gives earlier error (bioc 3.16 default version)
# install.packages("https://bioconductor.org/packages/3.16/bioc/src/contrib/treeio_1.22.0.tar.gz", repos=NULL, type="source")
## gives errors at the root step:
# ℹ invalid tbl_tree object. Missing column: parent.
# ℹ invalid tbl_tree object. Missing column: parent.
# ℹ invalid tbl_tree object. Missing column: parent.
# ℹ invalid tbl_tree object. Missing column: parent.
# ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
# ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

#### treeio_1.20.2 also gives earlier error. I think it worked before because I had an additional package loaded
install.packages("https://bioconductor.org/packages/3.15/bioc/src/contrib/treeio_1.20.2.tar.gz", repos=NULL, type="source")
library(treeio)
data(bird.orders, package="ape")
# add fake node labels
bird.orders$node.label <- paste("node",1:Nnode(bird.orders),sep="")
# convert to treedata
bird.orders_treedata <- as.treedata(bird.orders)

# confirm that there are meaningful tip and node labels
bird.orders_treedata@phylo$tip.label
bird.orders_treedata@phylo$node.label

##### reroot
bird.orders_treedata_rerooted <- root(bird.orders_treedata, "Galliformes")

# we no longer have those meaningful tip and node labels
bird.orders_treedata_rerooted@phylo$tip.label

bird.orders_treedata_rerooted@phylo$node.label



### here's the full set of packages I loaded before:
# library(Biostrings)
# library(ape)
# library(phylobase)
# # library(phangorn)
# library(ggtree)
# library(here)
# library(tidyverse)
# library(janitor)
# library(kableExtra)
# library(patchwork)
# library(ggmsa)  ## xx temp while I am having firewall trouble.  I think it's easy for me to isntall this on my Mac but not on the hutch server Rstudio.  But right now my mac isn't mounting the drive so I can't access the files
# library(tidytree)





### a minimal set of packages that I think could be involved. I WASN'T loading treeio directly before - one of the others must have been calling it.
# library(ape)
# library(phylobase)
# library(ggtree)
# library(tidytree)

### try various combinations:
library(treeio)
data(bird.orders, package="ape")
# add fake node labels
bird.orders$node.label <- paste("node",1:Nnode(bird.orders),sep="")
# convert to treedata
bird.orders_treedata <- as.treedata(bird.orders)

##### reroot
bird.orders_treedata_rerooted <- root(bird.orders_treedata, "Galliformes")
# ℹ invalid tbl_tree object. Missing column: parent.
# ℹ invalid tbl_tree object. Missing column: parent.
# ℹ invalid tbl_tree object. Missing column: parent.
# ℹ invalid tbl_tree object. Missing column: parent.
# ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
# ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

## restart R:
library(ape)
library(phylobase)
library(ggtree)
library(tidytree)
library(treeio)
data(bird.orders, package="ape")
# add fake node labels
bird.orders$node.label <- paste("node",1:Nnode(bird.orders),sep="")
# convert to treedata
bird.orders_treedata <- as.treedata(bird.orders)
##### reroot
bird.orders_treedata_rerooted <- root(bird.orders_treedata, "Galliformes")
bird.orders_treedata_rerooted 

## if only treeio is loaded:
# ℹ invalid tbl_tree object. Missing column: parent.
# ℹ invalid tbl_tree object. Missing column: parent.
# ℹ invalid tbl_tree object. Missing column: parent.
# ℹ invalid tbl_tree object. Missing column: parent.
# ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
# ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.

## if only tidytree is loaded: 
# Error in UseMethod("as.treedata") : 
#     no applicable method for 'as.treedata' applied to an object of class "phylo"

## if ggtree, tidytree and treeio are loaded, the root doesn't give error, but does lose the labels
## if ape, phylobase, ggtree, tidytree and treeio are loaded, the root doesn't give error, but does lose the labels


### xxx troubleshoot by comparing laptop to desktop computer?

