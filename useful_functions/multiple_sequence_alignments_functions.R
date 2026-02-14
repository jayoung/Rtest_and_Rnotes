### source this file as follows:
# malik_h_dir <- "/fh/fast/malik_h/"
# if (Sys.info()[["sysname"]]=="Darwin") { malik_h_dir <- "/Volumes/malik_h/" }
# source( paste0(malik_h_dir, "user/jayoung/git_more_repos/Rtest_and_Rnotes/useful_functions/multiple_sequence_alignments_functions.R") )

### degapAln - wrapper around maskGaps. Degaps alignment by columns
degapAln <- function(myAln, fractionOfSeqsWithGap=1) {
    maskedAln <- NULL
    if(class(myAln)=="AAStringSet") {
        maskedAln <- myAln |> 
            AAMultipleAlignment() |> 
            maskGaps(min.fraction=fractionOfSeqsWithGap, 
                     min.block.width=1) |> 
            AAStringSet()
    }
    if(class(myAln)=="DNAStringSet") {
        maskedAln <- myAln |> 
            DNAMultipleAlignment() |> 
            maskGaps(min.fraction=fractionOfSeqsWithGap, 
                     min.block.width=1) |> 
            DNAStringSet()
    }
    ## I don't think there is a BMultipleAlignment class!
    if(class(myAln)=="BStringSet") {
        maskedAln <- myAln |> 
            BMultipleAlignment() |> 
            maskGaps(min.fraction=fractionOfSeqsWithGap, 
                     min.block.width=1) |> 
            BStringSet()
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
                            frameshiftTranslatesTo="X", 
                            unknownCodonTranslatesTo="-", 
                            quiet=FALSE) {
    
    ## make new genetic code that includes codon="---" for which translation = "-"
    gapCodon <- "-"
    names(gapCodon) <- "---"
    my_GENETIC_CODE <- c(GENETIC_CODE, gapCodon)
    
    ## translate the codons
    pep <- my_GENETIC_CODE[myCodons]
    
    ## check for NA, from codons that were not possible to translate, e.g. frameshift codons and codons containing N
    if (sum(is.na(pep))>0) {
        ## maybe emit warnings
        if(!quiet) {
            myError <- paste0(
                "\nWarning - there were codons I could not translate. Using this character: ", 
                unknownCodonTranslatesTo, "\n")
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
        ## detect frameshift codons. The contain gap and non-gap characters:
        isFrameshift <- grepl("-", myCodons) & (nchar(gsub("-", "", myCodons)) > 0)
        pep[which(isFrameshift)] <- frameshiftTranslatesTo
        
        ## replace remaining NA with unknownCodonTranslatesTo
        pep[ which(is.na(pep)) ] <- unknownCodonTranslatesTo
    }
    
    ## prep for output
    pep <- paste(pep, collapse="")
    return(pep)
}

##### translateGappedAln - wrap the getCodons and translateCodons functions together into one:
# `frameshiftTranslatesTo="X"`    (e.g. for codons containing gap and DNA base), 
# `unknownCodonTranslatesTo="-"`  (e.g. for codons containing ambiguities) 
# `quiet=FALSE` : will emit warnings
translateGappedAln <- function(myAln, 
                               frameshiftTranslatesTo="X", 
                               unknownCodonTranslatesTo="-", 
                               quiet=FALSE) {
    ### some upfront checks
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
    
    ### maybe warn about codons containing ambiguities
    if(!quiet) {
        ## figure out warning about ambiguous codons
        ambiguityNucsPattern <- paste(setdiff(names(IUPAC_CODE_MAP), DNA_BASES), 
                                      collapse ="|")
        temp <- lapply(names(myCodons), function(x) {
            theseCodons <- myCodons[[x]]
            hasAmbig <- grepl(ambiguityNucsPattern, theseCodons)
            if(sum(hasAmbig)==0) { return(NULL) }
            ambigs <- theseCodons[which(hasAmbig)]
            ambigPositions <- which(hasAmbig)
            ambigMessage <- paste(ambigPositions, ambigs, sep="_")
            ambigMessage <- paste(ambigMessage, collapse=",")
            ambigMessage <- paste("\n\nWARNING - seq ",x," contains ambiguous codons:",ambigMessage,"\n\n", sep="")
            warning(ambigMessage)
        })
    }
    
    ### translate those codons
    myAAaln <- lapply(names(myCodons), function(x) {
        translateCodons(myCodons=myCodons[[x]], 
                        seqname=x,
                        frameshiftTranslatesTo=frameshiftTranslatesTo,
                        unknownCodonTranslatesTo=unknownCodonTranslatesTo, 
                        quiet=quiet)
    }) |> 
        set_names(names(myCodons)) |> 
        unlist() |> 
        AAStringSet()
    
    return(myAAaln)
}

########### alnToTibble - converts an alignment into a tibble, one row per alignment position, one column per aligned seq
## if we supply refSeqName, it'll also get the ungapped position in that refseq
alnToTibble <- function(aln, refSeqName=NULL) {
    names(aln) <- sapply(strsplit(names(aln), " "), "[[", 1)
    output <- as.matrix(aln) |> 
        t() |> 
        as_tibble() 
    output$aln_pos <- 1:width(aln)[1]
    output <- output |> 
        relocate(aln_pos)
    if(!is.null(refSeqName)) {
        if(!refSeqName %in% colnames(output)) {
            stop("\n\nERROR - the refSeqName you supplied is not in the alignment")
        }
        ref_seq <- output[,refSeqName] |> 
            deframe()
        output$ref_pos <- cumsum(ref_seq != "-")
        output[which(ref_seq=="-"),"ref_pos"] <- NA
        output <- output |> 
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
    }) |> 
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
    querySeqName <- unique(seqnames(queryCoords_gr)) |> 
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
    queryCoords_tbl <- queryCoords_gr |> 
        as.data.frame() |> 
        as_tibble()
    lookup <- lookup_tbl |> 
        select(all_of(c(querySeqName, targetSeqName)))
    colnames(lookup) <- c("query","target")
    
    queryCoords_tbl <- queryCoords_tbl |> 
        ## convert start
        left_join (lookup, by=c("start"="query")) |> 
        dplyr::rename(target_start=target) |> 
        ## convert end
        left_join (lookup, by=c("end"="query")) |> 
        dplyr::rename(target_end=target)
    new_gr <- queryCoords_tbl |> 
        select(start=target_start, end=target_end, strand=strand,
               orig_seqnames=seqnames, orig_start=start, orig_end=end)  |> 
        mutate(orig_width = orig_end + 1 - orig_start) |> 
        mutate(seqnames=targetSeqName) |> 
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
    querySeqNames <- unique(seqnames(queryCoords_gr)) |> 
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
    }) |> 
        unlist(use.names = FALSE)
    return(new_gr)
}



##### convertCoordsToOrigCoords takes output of convertCoords, after we're manipulated in in various ways, and returns a GRanges with the original coords
convertCoordsToOrigCoords <- function(gr) {
    gr |> 
        as.data.frame() |> 
        as_tibble() |> 
        select(-seqnames, -start, -end, -width, -orig_width)  |> 
        dplyr::rename_with( ~ str_remove_all(.x, "orig_")) |> 
        GRanges()
}


##### wrote this when working on Priya Shah YFV project (and Jyoti Batra bat genome comparisons)
### xxx I think the built-in pwalign::pid function gives the same answer as my pid_excl_gap column, but I want to check into that more.
### removeChars_forAlnLenCalc - sometimes I want to ignore positions containing certain characters (in addition to gaps) when I calculate degapped aln length
simple_percent_identity_pairwise <- function(twoSeqs_stringSet, 
                                             removeChars_forAlnLenCalc=NULL) {
    if(length(twoSeqs_stringSet) != 2) {
        stop("\n\nERROR - input does has the wrong number of sequences - function designed to work on just two sequences\n\n")
    }
    if(width(twoSeqs_stringSet)[1] !=  width(twoSeqs_stringSet)[2]) {
        stop("\n\nERROR - seqs don't have the same length as each other\n\n")
    }
    output <- tibble(seq1=names(twoSeqs_stringSet)[1],
                     seq2=names(twoSeqs_stringSet)[2] )
    
    twoSeqs <- as.matrix(twoSeqs_stringSet) |> 
        t() |> 
        as_tibble()
    colnames(twoSeqs) <- c("seq1", "seq2")
    
    twoSeqs_degapFullGap <- twoSeqs |> 
        filter(!(seq1=="-" & seq2=="-"))
    output$aln_len <- nrow(twoSeqs_degapFullGap)
    
    output$num_identical <- sum(twoSeqs_degapFullGap$seq1 == twoSeqs_degapFullGap$seq2)
    
    twoSeqs_degapAnyGap <- twoSeqs_degapFullGap |> 
        filter(seq1 != "-" & seq2 != "-")
    if(!is.null(removeChars_forAlnLenCalc)) {
        for (my_char in removeChars_forAlnLenCalc) {
            twoSeqs_degapAnyGap <- twoSeqs_degapFullGap |> 
                filter(seq1 != my_char & seq2 != my_char)
        }
    }
    
    output$aln_len_nogaps <- nrow(twoSeqs_degapAnyGap)
    
    output <- output |> 
        mutate(pid_incl_gap = 100*num_identical / aln_len) |> 
        mutate(pid_excl_gap = 100*num_identical / aln_len_nogaps)
    return(output)
}

simple_percent_identity_multiple <- function(one_aln) {
    if(length(unique(width(one_aln)))>1) {
        stop("\n\nERROR - these are not aligned seqs - they have different lengths: ", 
             paste(width(one_aln), collapse=","),
             "\n\n")
    }
    lapply(1:(length(one_aln)-1), function(i) {
        lapply( (i+1):length(one_aln), function(j) {
            simple_percent_identity_pairwise(one_aln[c(i,j)])
        }) |> 
            bind_rows()
    }) |> 
        bind_rows()
}

### xxx replace these examples with reproducible input
# tempAln <- alns[["pep"]][["janet"]][["Capsid"]]  # [9:11] #|> 
# narrow(start=1, end=20)
# tempAln |> 
#     simple_percent_identity_multiple()

# tempAln

# simple_percent_identity_multiple(alns[["pep"]][["tamanash"]][[1]])


####### alignedStringDistMatrix is a fiunction from Herve that is similar to stringDist but:
# (a) does not realign the seqs
# (b) allows supplying a matrix and an indel weight

# see https://github.com/Bioconductor/pwalign/issues/15

.normarg_weightmat <- function(weightmat, indel.weight=1L)
{
    if (!(is.matrix(weightmat) && is.numeric(weightmat)))
        stop(wmsg("'weightmat' must be NULL or a numeric matrix"))
    if (nrow(weightmat) != ncol(weightmat))
        stop(wmsg("'weightmat' must be a square matrix"))
    rownms <- rownames(weightmat)
    colnms <- colnames(weightmat)
    if (is.null(rownms) || is.null(colnms))
        stop(wmsg("'weightmat' must have rownames and colnames"))
    if (!(all(nchar(rownms) == 1L) && all(nchar(colnms) == 1L)))
        stop(wmsg("the rownames and colnames on 'weightmat' ",
                  "must be single letters"))
    if (!identical(rownms, colnms))
        stop(wmsg("the rownames and colnames on 'weightmat' ",
                  "must be the same"))
    if (!identical(t(weightmat), weightmat))
        stop(wmsg("'weightmat' must be symmetrical"))
    if ("-" %in% rownms) {
        if (!identical(indel.weight, 1L))
            warning("'indel.weight' is ignored when 'weightmat' has ",
                    "a row and a column for \"-\"", immediate.=TRUE)
        return(weightmat)
    }
    if (!isSingleNumber(indel.weight))
        stop(wmsg("'indel.weight' must be a single number"))
    ## Add "-" row and column to 'weightmat'.
    weightmat <- rbind(weightmat, `-`=indel.weight)
    weightmat <- cbind(weightmat, `-`=indel.weight)
    weightmat["-", "-"] <- 0L
    weightmat
}

### Compute the distance matrix of a set of aligned strings using the
### Hamming distance or weighted Hamming distance.
###
### - 'x' must be a character vector (or XStringSet derivative like
###   DNAStringSet or AAStringSet) where all strings have the same length.
###
### - 'weightmat' must be a square symmetric matrix that contains the
###   weights of the matches and mismatches between 2 given letters. Similar
###   to a scoring matrix like BLOSUM62 except that weights are conceptually
###   the opposite of scores. This means that when using BLOSUM62 with
###   alignedStringDistMatrix(), remember to pass -BLOSUM62.
###
### - 'indel.weight' is the weight of a mismatch between "-" and any other
###   letter. By default the weight of a match between "-" and itself is
###   considered to be 0. Note that 'indel.weight' will be ignored
###   if 'weightmat' has a row and a column for "-".
###
### Returns a square symmetric matrix with length(x) rows and columns.

alignedStringDistMatrix <- function(x, weightmat=NULL, indel.weight=1L)
{
    if (!is(x, "BStringSet"))
        x <- BStringSet(x)
    if (!isConstant(width(x)))
        stop(wmsg("all the strings in 'x' must have the same length"))
    if (is.null(weightmat)) {
        ## Computes the weights between 2 vectors of letters of the same
        ## length. Weight is 0 for a match, 1 for a mismatch that doesn't
        ## involve "-", and 'indel.weight' for a mismatch that involves "-".
        ## Returns a numeric vector parallel to 'vol1' and 'vol2'.
        get_weigths <- function(vol1, vol2) {
            not_equal <- vol1 != vol2
            is_minus <- vol1 == "-" | vol2 == "-"
            (not_equal & !is_minus) + (not_equal & is_minus) * indel.weight
        }
    } else {
        weightmat <- .normarg_weightmat(weightmat, indel.weight)
        if (!all(uniqueLetters(x) %in% rownames(weightmat)))
            stop(wmsg("'x' contains letters that are not in the dimnames ",
                      "of 'weightmat'"))
        ## Computes the weights between 2 vectors of letters of the
        ## same length. The weight between 2 given letters is obtained
        ## from 'weightmat'.
        ## Returns a numeric vector parallel to 'vol1' and 'vol2'.
        get_weigths <- function(vol1, vol2) weightmat[cbind(vol1, vol2)]
    }
    ans <- matrix(NA_integer_, nrow=length(x), ncol=length(x))
    x_names <- names(x)
    if (!is.null(x_names))
        dimnames(ans) <- list(x_names, x_names)
    ## Turn 'x' into a list of vectors of letters.
    vol_list <- lapply(x, function(xx) rawToChar(as.raw(xx), multiple=TRUE))
    for (i in seq_along(vol_list)) {
        vol1 <- vol_list[[i]]
        for (j in seq_len(i)) {
            vol2 <- vol_list[[j]]
            ans[i, j] <- ans[j, i] <- sum(get_weigths(vol1, vol2))
        }
    }
    ans
}



######## these functions are from jensenShannonDistance.Rmd (and also defined there, so I don't want to change one without the other. I'm sure there's a better way to handle that)

# a tiny function that makes sure all seqs in an alignment are the same length as each other
checkAlnLengths <- function(aln) {
    if(length(unique(width(aln))) != 1) {
        stop("\n\nERROR - you supplied a ragged alignment (seqs not all the same length)\n\n")
    } else {
        return(TRUE)
    } 
}

## define the letters we want to count
# AA_STANDARD is defined in the Biostrings package and includes the usual 20 amino acids. I want to add the gap character ("-")
myAAtoTabulate <- c("-", AA_STANDARD)

## define the function
getAlnCounts <- function(aln, letters=myAAtoTabulate, as.prob=FALSE) {
    # check for ragged alns (seqs not all the same length)
    checkAlnLengths(aln)
    
    ## for each sequence, get matrix of 0 and 1 representing the letter at each position. Returns a list of matrices, one for each input seq, with num_rows= num aln positions, and num_columns=21 (- plus each AA)  
    countsEachSeq <- lapply(1:length(aln), function(i) {
        letterFrequencyInSlidingView(aln[[i]], view.width = 1, letters=letters)
    })
    
    # if there were letters in the alignment that are not accounted for in the letters argument, the totals won't be correct.
    expectedTotals <- width(aln)[1]
    totalCountsEachSeq <- sapply(countsEachSeq, sum)
    if ( sum(totalCountsEachSeq != expectedTotals) > 0) {
        ## generate an informative error message:
        problem_seqs <- which(totalCountsEachSeq != expectedTotals)
        problem_seq_letters <- aln[problem_seqs] |> 
            as.character() |> 
            strsplit(split="")
        problem_seq_letters <- lapply(problem_seq_letters, function(x) {
            setdiff(x, myAAtoTabulate) |> 
                unique() |> 
                paste0()
        }) |> 
            paste(collapse=",")
        err_msg <- paste0("\n\nERROR - the total counts didn't add up correctly.\n",
                          "These sequences contain unexpected letters: " , 
                          paste(problem_seqs, collapse=","),
                          "\nAnd those letters are: ",
                          problem_seq_letters,
                          "\n\n")
        stop(err_msg)
    }
    
    # get total counts by position - the Reduce function takes a list object and uses the specified function on all the elements
    countTotals <- Reduce("+", countsEachSeq)
    
    # transpose so columns are positions and rows are each letter type
    countTotals <- t(countTotals)
    
    # perhaps get frequencies not counts
    if(as.prob) {
        freqs <- countTotals / colSums(countTotals)
        return(freqs)
    } else {
        return(countTotals)
    }
}

