### source this file as follows:
# malik_h_dir <- "/fh/fast/malik_h/"
# if (Sys.info()[["sysname"]]=="Darwin") { malik_h_dir <- "/Volumes/malik_h/" }
# source( paste0(malik_h_dir, "user/jayoung/git_more_repos/Rtest_and_Rnotes/useful_functions/multiple_sequence_alignments_functions.R") )

### degapAln - wrapper around maskGaps. Degaps alignment by columns
degapAln <- function(myAln, fractionOfSeqsWithGap=1) {
    maskedAln <- NULL
    if(class(myAln)=="AAStringSet") {
        maskedAln <- myAln %>% 
            AAMultipleAlignment() %>% 
            maskGaps(min.fraction=fractionOfSeqsWithGap, 
                     min.block.width=1) %>% 
            AAStringSet()
    }
    if(class(myAln)=="DNAStringSet") {
        maskedAln <- myAln %>% 
            DNAMultipleAlignment() %>% 
            maskGaps(min.fraction=fractionOfSeqsWithGap, 
                     min.block.width=1) %>% 
            DNAStringSet()
    }
    if(is.null(maskedAln)) {
        stop("\n\nERROR - input alignment is an object of a class that the function is not set up to handle\n\n")
    }
    return(maskedAln)
}

########  codon_range_to_nuc_range - given a codon position range, return the equivalent nucleotide position range to take. This function is for simple pairs of coordinates. 
codon_range_to_nuc_range <- function(codon_start, codon_end=codon_start) {
    if(codon_start>codon_end) {
        stop("\n\nERROR - start is after end\n\n")
    }
    nuc_start <- (codon_start-1)*3 + 1
    nuc_end <- (codon_end-1)*3 + 3
    return(c(nuc_start, nuc_end))
}
## a version for GRanges
codon_range_to_nuc_range_gr <- function(gr) {
    if(class(gr) != "GRanges") {
        stop("\n\nERROR - this function is designed for GRanges objects\n\n")
    }
    nuc_starts <- (start(gr)-1)*3 + 1
    nuc_ends <- (end(gr)-1)*3 + 3
    new_gr <- gr
    ranges(new_gr) <- IRanges(start=nuc_starts, end=nuc_ends)
    return(new_gr)
}


########  nuc_range_to_codon_range - given a nucleotide position range, return the equivalent amino acid position range to take.
nuc_range_to_codon_range <- function(nuc_start, nuc_end=nuc_start) {
    if(nuc_start>nuc_end) {
        stop("\n\nERROR - start is after end\n\n")
    }
    if((nuc_end+1-nuc_start) %% 3 != 0) {
        stop("\n\nERROR - range is not an even multiple of 3\n\n")
    }
    codon_start <- ceiling(nuc_start/3)
    codon_end <- ceiling(nuc_end/3)
    return(c(codon_start, codon_end))
}

## a version for GRanges
nuc_range_to_codon_range_gr <- function(gr) {
    if(class(gr) != "GRanges") {
        stop("\n\nERROR - this function is designed for GRanges objects\n\n")
    }
    testWidths <- width(gr) %% 3
    if (sum(testWidths != 0) > 0) {
        problemRanges <- which(testWidths != 0)
        problemRanges <- gr[problemRanges]
        cat("\n### problems!\n")
        print(problemRanges)
        stop("\n\nERROR - there are ranges that are not an even multiple of 3 wide\n\n")
    }
    codon_start <- ceiling(start(gr)/3)
    codon_end <- ceiling(end(gr)/3)
    
    new_gr <- gr
    ranges(new_gr) <- IRanges(start=codon_start, end=codon_end)
    return(new_gr)
}


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
                            seqname=NULL,
                            unknownCodonTranslatesTo="-", 
                            # ambiguousNucHandling="warn",
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
            myError <- paste0("\nWarning - there were codons I could not translate. Using this character: ", unknownCodonTranslatesTo, "\n")
            if(!is.null(seqname)) {
                myError <- gsub("Warning - there were",
                                paste0("Warning in seq ",seqname," - there were"),
                                myError)
            }
            message(myError)
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

##### translateGappedAln - wrap the getCodons and translateCodons functions together into one:
# ambiguousNucHandling - this is for codons containing N/Y/etc
#    none:   replaces with unknownCodonTranslatesTo
#    warn:   replaces with unknownCodonTranslatesTo and emits warning
#    gap:    replaces with a gap codon
#    stop:   throws an error
translateGappedAln <- function(myAln, 
                               unknownCodonTranslatesTo="-", 
                               ambiguousNucHandling="warn",
                               quiet=FALSE) {
    ### some upfront checks
    if(!ambiguousNucHandling %in% c("none", "warn", "stop", "gap")) {
        stop("\n\nthe ambiguousNucHandling option must be either none, warn or stop\n\n")
    }
    if(length(names(myAln)) != length(unique(names(myAln)))) {
        duplicatedNames <- names(myAln)
        duplicatedNames <- duplicatedNames[which(duplicated(duplicatedNames))]
        duplicatedNames <- unique(duplicatedNames)
        duplicatedNames <- paste(duplicatedNames, collapse=",")
        stop("\n\nERROR - the alignment does not have unique names: the duplicates are: ",
             duplicatedNames,"\n\n")
    }
    
    ### split seqs into codons (list of character vectors)
    myCodons <- getCodons(myAln)

    ### deal with codons containing ambiguities
    if(ambiguousNucHandling != "none") {
        ambiguityNucsPattern <- paste(setdiff(names(IUPAC_CODE_MAP), DNA_BASES), 
                                      collapse ="|")
        myCodons <- lapply(names(myCodons), function(x) {
            theseCodons <- myCodons[[x]]
            hasAmbig <- grepl(ambiguityNucsPattern, theseCodons)
            if(sum(hasAmbig)==0) { return(theseCodons) }
            ## figure out warning about amiguous codons
            ambigs <- theseCodons[which(hasAmbig)]
            ambigPositions <- which(hasAmbig)
            ambigMessage <- paste(ambigPositions, ambigs, sep="_")
            ambigMessage <- paste(ambigMessage, collapse=",")
            ambigMessage <- paste("\n\nWARNING - seq ",x," contains ambiguous codons:",ambigMessage,"\n\n", sep="")
            
            if(ambiguousNucHandling == "warn") {
                if(!quiet) {warning(ambigMessage)}
            }
            if(ambiguousNucHandling == "stop") {
                ambigMessage <- gsub("\n\n$",
                                     " (and maybe other seqs, don't know)\n\n",
                                     ambigMessage)
                stop(ambigMessage) 
            }
            if(ambiguousNucHandling == "gap") {
                theseCodons[which(hasAmbig)] <- "---"
            }
            return(theseCodons)
        }) %>% 
            set_names(names(myCodons))
    }

    ### translate those codons
    myAAaln <- lapply(names(myCodons), function(x) {
        translateCodons(myCodons=myCodons[[x]], 
                        seqname=x,
                        unknownCodonTranslatesTo=unknownCodonTranslatesTo, 
                        quiet=quiet)
    }) %>% 
        set_names(names(myCodons)) %>% 
        unlist() %>% 
        AAStringSet()
    
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



