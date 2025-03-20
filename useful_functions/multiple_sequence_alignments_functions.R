### source this file as follows:
# malik_h_dir <- "/fh/fast/malik_h/"
# if (Sys.info()[["sysname"]]=="Darwin") { malik_h_dir <- "/Volumes/malik_h/" }
# source( paste0(malik_h_dir, "user/jayoung/git_more_repos/Rtest_and_Rnotes/useful_functions/multiple_sequence_alignments_functions.R") )



####### translateGappedAln:  translate() does not deal with gap characters yet (should later - see https://github.com/Bioconductor/Biostrings/issues/30). So I made my own translateGappedAln function.

## getCodons - splits alignment into constituent codons. Utility function for translateGappedAln
getCodons <- function(myAln) {
    seqs <- as.character(myAln)
    len <- width(myAln)[1]
    starts <- seq(from=1, to=len, by=3)
    ends <- starts + 2
    myViews <- lapply(myAln, function(x) { 
        Views(x, starts, ends)
    })
    myCodons <- lapply(myViews, function(x) {
        as.character(DNAStringSet(x))
    })
    myCodons
}

## translateCodons - takes a character vector of codons as input, outputs the corresponding amino acids. Utility function for translateGappedAln
translateCodons <- function(myCodons, 
                            unknownCodonTranslatesTo="-", 
                            quiet=FALSE) {
    ## make new genetic code
    gapCodon <- "-"
    names(gapCodon) <- "---"
    my_GENETIC_CODE <- c(GENETIC_CODE, gapCodon)
    
    ## translate the codons
    pep <- my_GENETIC_CODE[myCodons]
    
    ## check for codons that were not possible to translate, e.g. frameshift codons
    if (sum(is.na(pep))>0) {
        if(!quiet) {
            message("\nWarning - there were codons I could not translate. Using this character: ", unknownCodonTranslatesTo, "\n")
            unknownCodons <- unique(myCodons[ which(is.na(pep))])
            message("The codons in question were: ",
                    paste(unknownCodons, collapse=","),
                    "\n\n")
        }
        pep[ which(is.na(pep)) ] <- unknownCodonTranslatesTo
    }
    
    ## prep for output
    pep <- paste(pep, collapse="")
    return(pep)
}

## wrap the getCodons and translateCodons functions together into one:
translateGappedAln <- function(myAln, 
                               unknownCodonTranslatesTo="-", 
                               quiet=FALSE) {
    myCodons <- getCodons(myAln)
    myAAaln <- AAStringSet(unlist(lapply(
        myCodons, 
        translateCodons, 
        unknownCodonTranslatesTo=unknownCodonTranslatesTo, 
        quiet=quiet)))
    return(myAAaln)
}


########### alnToTibble - converts an alignment into a tibble, one row per alignment position, one column per aligned seq
## if we supply refSeqName, it'll also get the ungapped position in that refseq
alnToTibble <- function(aln, refSeqName=NULL) {
    names(aln) <- sapply(strsplit(names(aln), " "), "[[", 1)
    output <- as.matrix(aln) %>% 
        t() %>% 
        as_tibble() 
    output$aln_pos <- 1:width(aln)[1]
    output <- output %>% 
        relocate(aln_pos)
    if(!is.null(refSeqName)) {
        if(!refSeqName %in% colnames(output)) {
            stop("\n\nERROR - the refSeqName you supplied is not in the alignment")
        }
        ref_seq <- output[,refSeqName] %>% 
            deframe()
        output$ref_pos <- cumsum(ref_seq != "-")
        output[which(ref_seq=="-"),"ref_pos"] <- NA
        output <- output %>% 
            relocate(ref_pos, .after=aln_pos)
    }
    return(output)
}


########### getAlnPosLookupTable - uses an alignment to get a lookup table of alignment position versus ungapped position in each sequence. Gap positions get an NA
getAlnPosLookupTable <- function(myAln) {
    if(length(unique(width(myAln)))>1) {
        stop("\n\nERROR - the seqs in the alignment are not all the same length\n\n")
    }
    output <- tibble(aln_pos = 1:width(myAln)[1])
    myAln_chars <- strsplit(as.character(myAln),"")
    each_seq_pos <- sapply(myAln_chars, function(eachSeq_chars) {
        countNonGap <- cumsum(eachSeq_chars != "-")
        countNonGap[which(eachSeq_chars=="-")] <- NA
        return(countNonGap)
    }) %>% 
        as_tibble()
    output <- bind_cols(output, each_seq_pos)
    return(output)
}




##### convertCoordsOneSeq - a function, given a lookup table produced by getAlnPosLookupTable, and a GRanges of positions in one aligned sequence, to make a GRanges of the equivalent positions in another seq, or in the alignment 
## queryCoords_gr is a GRanges containing only one unique seqnames
convertCoordsOneSeq <- function(queryCoords_gr, 
                                lookup_tbl, 
                                targetSeqName="aln_pos", 
                                keep_mcols=TRUE) {
    querySeqName <- unique(seqnames(queryCoords_gr)) %>% 
        as.character()
    if(length(querySeqName) > 1) {
        stop("\n\nERROR in convertCoordsOneSeq function - the GRanges objects contains multiple different seqnames\n\n")
    }
    if(!querySeqName %in% names(lookup_tbl)) {
        stop("\n\nERROR in convertCoordsOneSeq function - the query seq ", querySeqName, " is not in the lookup table colnames\n\n")
    }
    if(!targetSeqName %in% colnames(lookup_tbl)) {
        stop("\n\nERROR in convertCoordsOneSeq function - the target seq ", targetSeqName, " is not in the lookup table colnames\n\n")
    }
    queryCoords_tbl <- queryCoords_gr %>% 
        as.data.frame() %>% 
        as_tibble()
    lookup <- lookup_tbl %>% 
        select(all_of(c(querySeqName, targetSeqName)))
    colnames(lookup) <- c("query","target")
    
    queryCoords_tbl <- queryCoords_tbl %>% 
        ## convert start
        left_join (lookup, by=c("start"="query")) %>% 
        dplyr::rename(target_start=target) %>% 
        ## convert end
        left_join (lookup, by=c("end"="query")) %>% 
        dplyr::rename(target_end=target)
    new_gr <- queryCoords_tbl %>% 
        select(start=target_start, end=target_end, strand=strand,
               orig_seqnames=seqnames, orig_start=start, orig_end=end)  %>% 
        mutate(orig_width = orig_end + 1 - orig_start) %>% 
        mutate(seqnames=targetSeqName) %>% 
        GRanges()
    if(keep_mcols) { mcols(new_gr) <- cbind(mcols(new_gr),
                                            mcols(queryCoords_gr)) }
    return(new_gr)
}
# convertCoordsOneSeq(temp_orfs, genomes_aln_02_lookup)


##### convertCoords is a wrapper function to run convertCoordsOneSeq on GRanges containing >1 different query sequence 
convertCoords <- function(queryCoords_gr, 
                          lookup_tbl, 
                          targetSeqName="aln_pos", 
                          keep_mcols=TRUE) {
    querySeqNames <- unique(seqnames(queryCoords_gr)) %>% 
        as.character()
    missingQueryNames <- setdiff(querySeqNames, colnames(lookup_tbl))
    if(length(missingQueryNames) > 0) {
        stop("\n\nERROR in convertCoords function - the GRanges objects contains seqnames that are not in the lookup table: ",missingQueryNames, "\n\n")
    }
    if(!targetSeqName %in% colnames(lookup_tbl)) {
        stop("\n\nERROR in convertCoordsOneSeq function - the target seq ", targetSeqName, " is not in the lookup table colnames\n\n")
    }
    queryCoords_gr_by_seq <- split(queryCoords_gr, 
                                   as.character(seqnames(queryCoords_gr)))
    new_gr <- endoapply(queryCoords_gr_by_seq, function(x) {
        convertCoordsOneSeq(x, 
                            lookup_tbl=lookup_tbl, 
                            targetSeqName=targetSeqName, 
                            keep_mcols=keep_mcols)
    }) %>% 
        unlist(use.names = FALSE)
    return(new_gr)
}



##### convertCoordsToOrigCoords takes output of convertCoords, after we're manipulated in in various ways, and returns a GRanges with the original coords
convertCoordsToOrigCoords <- function(gr) {
    gr %>% 
        as.data.frame() %>% 
        as_tibble() %>% 
        select(-seqnames, -start, -end, -width, -orig_width)  %>% 
        dplyr::rename_with( ~ str_remove_all(.x, "orig_")) %>% 
        GRanges()
}


