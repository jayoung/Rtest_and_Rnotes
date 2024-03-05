trees_usingPublicTrees
================
Janet Young

2024-02-09

``` r
## get example species names I want a tree for. This is long-winded but it should be a working example:
## running from rhino
# RTL3_alnFile <- "/fh/fast/malik_h/user/jayoung/miscMalikLab/domesticated_capsid/alignments/all_alignment_dirs/RTL3/RTL3_aln2022summer_v28.fa"
## running from mac:
RTL3_alnFile <- "/Volumes/malik_h/user/jayoung/miscMalikLab/domesticated_capsid/alignments/all_alignment_dirs/RTL3/RTL3_aln2022summer_v28.fa"

RTL3_aln <- readDNAStringSet(RTL3_alnFile)
names(RTL3_aln) <- sapply(strsplit(names(RTL3_aln), " "), "[[", 1)

## now make a seqinfo table - seqname, seqlength, these overlaps, num capsid/protease-overlapping stop-free regions
RTL3_alnInfoTable <- tibble(id=names(RTL3_aln))
## get species name from id
RTL3_alnInfoTable <- RTL3_alnInfoTable %>% 
    separate(id, into=c(NA,"species"), extra="drop", remove=FALSE) %>% 
    mutate(species=str_replace(species, "CapsidHMMhit", "")) %>% 
    mutate(species=str_replace(species, "ProteaseORF", "")) %>% 
    mutate(species=str_replace(species, "BetweenORFsRegion", "")) %>% 
    mutate(species=str_to_sentence(to_snake_case(species))) %>% 
    mutate(species=str_replace_all(species,"_"," "))
```

``` r
specNamesToGetTreeFor <- RTL3_alnInfoTable %>% 
    # filter(!is.na(name_in_fig)) %>% 
    select(species) %>% 
    deframe() %>% 
    unique()
```

``` r
## add commonNameLookup to RTL3_alnInfoTable
commonNameLookup_tbl <- commonNameLookup %>% 
    unlist() %>% 
    as_tibble(rownames="species")
commonNameLookup_tbl <- commonNameLookup_tbl %>% 
    dplyr::rename("common_name"="value")
RTL3_alnInfoTable <- left_join(RTL3_alnInfoTable, commonNameLookup_tbl, by="species") %>% 
    relocate(common_name, .after="species")
```

``` r
### try taxize package

## try getting a tree, small test:
specNamesToGetTreeFor_test <- sample(specNamesToGetTreeFor, 8)

## classification() returns a list of data.frames, one per species queried
# tested what happens when a query species is not found: we get NA (logical) back for that list element
specNamesToGetTreeFor_test_classification <- classification(c(specNamesToGetTreeFor_test, "fake species name", "fake species name2"),
                                                            db="ncbi")
# numSpeciesNotFound <- sum(sapply(specNamesToGetTreeFor_test_classification, function(x) {!class(x)=="data.frame"}))
# and class2tree returns a tree of the non-NA elements
specNamesToGetTreeFor_test_classification_tree <- class2tree(specNamesToGetTreeFor_test_classification)

## get full tree - NCBI classifications
specNames_classification_ncbi <- classification(specNamesToGetTreeFor, db="ncbi")
# check how many species not found:
# numSpeciesNotFound <- sum(sapply(specNames_classification, function(x) {!class(x)=="data.frame"}))
# which(sapply(specNames_classification, function(x) {!class(x)=="data.frame"}))
## using ncbi database, all are found. Tree is not always well-resolved
specNames_tree_ncbi <- class2tree(specNames_classification_ncbi)
plot(specNames_tree_ncbi)


## get full tree - encyclopedia of life classifications
temp2 <- get_eolid("Pinus contorta", ask=FALSE, rows=1) 
temp1 <- get_eolid("Pinus contorta") 
get_eolid(specNamesToGetTreeFor[1:2], ask=FALSE, rows=1) 
specNames_classification_eol <- classification(specNamesToGetTreeFor[1:2], db="eol", rows=1)
# found both but then
# Error: Bad Gateway (HTTP 502)
# and no result

## get full tree - Integrated Taxonomic Information Service (itis) classifications

# without the rows=3, it stalls on "Chinchilla lanigera". For "Chinchilla lanigera", rows=1 was not enough to get info. Turns out there are 2 tsn (like taxID?) for Chinchilla, first one invalid. Using rows>=2 allows the second (valid) ID to be used

specNames_classification_itis <- classification(specNamesToGetTreeFor, db="itis", accepted=TRUE)


## Ochotona curzoniae needs rows=1 (or NA)
classification("Ochotona curzoniae", db="itis", rows=1)
classification("Ochotona curzoniae", db="itis", rows=2)


## Chinchilla lanigera needs rows=2 (or more)
classification("Chinchilla lanigera", db="itis", rows=1)
classification("Chinchilla lanigera", db="itis", rows=2)


# aha! this works
classification("Chinchilla lanigera", db="itis", accepted=TRUE)

itis_acceptname("Chinchilla lanigera")

temp_itis <- classification("Ochotona curzoniae", db="itis", rows=10)



## check how many species not found:
# table(sapply(specNames_classification_itis, class))
# numSpeciesNotFound <- sum(sapply(specNames_classification_itis, function(x) {!class(x)=="data.frame"}))
# which(sapply(specNames_classification_itis, function(x) {!class(x)=="data.frame"}))

# Nannospalax galili Cricetulus griseus 
classification("Nannospalax galili", db="itis") # not found
classification("Cricetulus griseus", db="itis")

classification("Nannospalax galili", db="ncbi") 
classification("Cricetulus griseus", db="ncbi")
temp1 <- classification("Homo sapiens", db="itis", accepted=TRUE)
class(temp1[["Homo sapiens"]])

temp2 <- classification("Nannospalax galili", db="itis", accepted=TRUE)
class(temp2[["Nannospalax galili"]])

### use itis when possible, if not, use NCBI
# there's more advice here: https://docs.ropensci.org/taxize/articles/name_cleaning.html
getClassification <- function(oneSpecName) {
    thisClassification <- classification(oneSpecName, db="itis", accepted=TRUE)
    # return(thisClassification)
    if(class(thisClassification[[oneSpecName]]) != "data.frame") {
        cat ("   using NCBI for species ",oneSpecName,"\n")
        thisClassification <- classification(oneSpecName, db="ncbi", accepted=TRUE)
    }
    return(thisClassification)
}
specNames_classification_itisNCBI <- lapply(specNamesToGetTreeFor, getClassification)
names(specNames_classification_itisNCBI) <- specNamesToGetTreeFor
## check how many species not found - all are now found
# table(sapply(specNames_classification_itisNCBI, class))
# table(sapply(specNames_classification_itisNCBI, function(x) {class(x[[1]])}))
# numSpeciesNotFound <- sum(sapply(specNames_classification_itisNCBI, function(x) {!class(x)=="data.frame"}))
# which(sapply(specNames_classification_itisNCBI, function(x) {!class(x)=="data.frame"}))

# myDim <- t(sapply(specNames_classification_itisNCBI, function(x) {dim(x[[1]])}))

specNames_classification_itisNCBI_reformat <- lapply(specNames_classification_itisNCBI, "[[", 1)
attr(specNames_classification_itisNCBI_reformat,"class") <- "classification"
attr(specNames_classification_itisNCBI_reformat,"db") <- "ncbi"


## using ncbi database, all are found. Tree is not always well-resolved
specNames_tree_itisNCBI <- class2tree(specNames_classification_itisNCBI_reformat)
plot(specNames_tree_itisNCBI)
```