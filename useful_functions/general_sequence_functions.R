### source this file as follows:
# malik_h_dir <- "/fh/fast/malik_h/"
# if (Sys.info()[["sysname"]]=="Darwin") { malik_h_dir <- "/Volumes/malik_h/" }
# source( paste0(malik_h_dir, "user/jayoung/git_more_repos/Rtest_and_Rnotes/useful_functions/general_sequence_functions.R") )

###### read_seq_file_JY is a utility function to read a fasta file
## reads file, and parses out seq IDs and description lines
## returns a list of the sequences (as DNAStringSet etc) and a tibble of seq_ids and descriptions
## type can be DNA or AA or B (BStringSet)
read_seq_file_JY <- function(seq_file, 
                             type="DNA", 
                             strip_afterDot=FALSE,
                             nrec = -1L  # for testing purposes, read the first few records
) {
    if(!file.exists(seq_file)) {
        stop("\n\nERROR - seqfile does not exist: ",seq_file, "\n\n")
    }
    if (type=="DNA") {
        seqs <- seq_file |> readDNAStringSet(nrec=nrec)
    }
    if (type=="AA") {
        seqs <- seq_file |> readAAStringSet(nrec=nrec)
    }
    if (type=="B") {
        seqs <- seq_file |> readBStringSet(nrec=nrec)
    }
    if(!exists("seqs")) {
        stop("\n\nERROR - you specified an invalid type. Should be DNA, AA, or B\n\n")
    }
    headers_split <- names(seqs) |> strsplit(" ")
    info_tbl <- tibble(seq_id = sapply(headers_split, "[[", 1),
                       seq_len = width(seqs),
                       desc = sapply(headers_split, function(x) {
                           if(length(x)>1) {
                               return( paste(x[2:length(x)], collapse=" ")  )
                           } else { return(NA) }
                       }))
    if(strip_afterDot) {
        info_tbl <- info_tbl |> 
            mutate(seq_id_withDot = seq_id) |> 
            mutate(seq_id = str_remove_all(seq_id, "\\.\\d+$")) |> 
            relocate(seq_id_withDot, .after=seq_id)
    }
    if (length(info_tbl$seq_id) == length(unique(info_tbl$seq_id))) {
        names(seqs) <- info_tbl$seq_id
    } else {
        warning("\nWARNING - not stripping descriptions off the seqnames, because then the IDs will not be unique\n\n")
    }
    return(list(seqs=seqs, info=info_tbl))
}
# fer_ensembl_cds_file <- "/fh/fast/malik_h/grp/public_databases/Ensembl/release-115/Rhinolophus_ferrumequinum/Rhinolophus_ferrumequinum.mRhiFer1_v1.p.cds.all.fa"
# test <- read_seq_file_JY(fer_ensembl_cds_file, type="DNA", strip_afterDot=TRUE, nrec=10)
# test[["seqs"]] |> names()
# test[["info"]]

####### read_multiple_seq_files_JY - when we have >1 seqfile, uses read_seq_file_JY to read them, and combines the results, add a "file" column to the info table
read_multiple_seq_files_JY <- function(files, ...) {
    dat <- lapply(files, read_seq_file_JY, ...)
    files_without_dir <- strsplit(files, "/")
    files_without_dir <- sapply(files_without_dir, function(x) {
        x[length(x)]
    })
    names(dat) <- files_without_dir
    
    ## seqs
    seqs <- lapply(dat, "[[", "seqs")
    if(class(seqs[[1]]) == "DNAStringSet") {
        seqs <- seqs |> DNAStringSetList() |> unlist(use.names=FALSE)
    }
    if(class(seqs[[1]]) == "AAStringSet") {
        seqs <- seqs |> AAStringSetList() |> unlist(use.names=FALSE)
    }
    if(class(seqs[[1]]) == "BStringSet") {
        seqs <- seqs |> BStringSetList() |> unlist(use.names=FALSE)
    }
    ## info
    info <- lapply(names(dat), function(x) {
        dat[[x]][["info"]] |> 
            mutate(file=x) |>
            relocate(file)
    }) |> 
        bind_rows()
    output <- list(seqs=seqs, info=info)
    return(output)
}

#### countAmbiguitiesEachSeq - a small function that counts all non-ACGT seqs in each member of a DNAStringSet
countAmbiguitiesEachSeq <- function(dna_seqs) {
    alphabetFrequency(dna_seqs) |> 
        as_tibble() |> 
        select(-A, -C, -G, -T) |> 
        rowSums()
}

###### GCcontent - returns GC content values on a single sequence which is a Biostrings object
## returns num GCs / num ACGTs (i.e. ignores other nucleotides in numerator and denominator)
GCcontent <- function (myseqs) {
    if (class(myseqs)=="DNAString") { myseqs <- DNAStringSet(myseqs) }
    if (class(myseqs)!="DNAStringSet") { 
        stop("\n\nERROR - input needs to be either a DNAString or a DNAStringSet object\n\n")
    }
    
    nuc_counts <- alphabetFrequency(myseqs)[,c("A","C","G","T")]
    ## in the case where we only had one seq, nuc_counts is an integer vector, otherwise it's a matrix. 
    if(class(nuc_counts)[1]=="integer") {
        gc_totals <- sum(nuc_counts[c("C","G")])
        grand_totals <- sum(nuc_counts)
    } else {
        gc_totals <- rowSums(nuc_counts[,c("C","G")])
        grand_totals <- rowSums(nuc_counts)
    }
    gc_fraction <- gc_totals / grand_totals
    if(!is.null(names(myseqs))) {
        names(gc_fraction) <- names(myseqs)
    }
    return(gc_fraction)
}
## ## test code 
# myseqs <- c(seq1="AGTAGTGCATGTATGC",
#             seq2="GGTAGCTGGATGATCGGTA") |>
#     DNAStringSet()
# # alphabetFrequency(myseqs)
# # width(myseqs)
# 
# GCcontent(myseqs[1])
# GCcontent(myseqs[2])
# GCcontent(myseqs)
## correct answers: 0.4375, 0.526





####### readBlast reads a processed blast file (processed via blastparsenew.bioperl )
# I include numLinesToRead option so I can do smaller tests on large files
readBlast <- function(procBlastFile, numLinesToRead=Inf, quiet=FALSE) {
    if (!file.exists(procBlastFile)) {
        stop("\n\nERROR - input file doesn't exist: ",procBlastFile,"\n\n")
    }
    ## first, check for the case where there were absolutely no hits in blast
    first_line_fields <- procBlastFile |> 
        scan(what="character", sep="\n", nmax=1, quiet = TRUE) |> 
        strsplit("\t")
    first_line_fields <- first_line_fields[[1]]
    if (first_line_fields[2] == "no_hits_in_blast") {
        if(!quiet) {
            cat("    WARNING - absolutely no blast hits in file",
            procBlastFile,
            " - it won't appear in the output tibble\n\n")
        }
        return(NULL)
    }
    
    ## get header
    header <- read_delim(procBlastFile, n_max=3, show_col_types=FALSE)
    
    ## get the rest of the data
    dat <- read_delim(procBlastFile, skip=2, col_names=FALSE, 
                      show_col_types=FALSE, n_max=numLinesToRead)
    colnames(dat) <- colnames(header)
    dat <- dat %>% clean_names()
    ## deal with the case where there were no hits with good enough score to make it into the parsed blast file:
    if(nrow(dat)==0) {
        dat <- header |> 
            clean_names() |> 
            filter(query != "-----")
        ## when there's no actual data, read_delim doesn't understand what class to make the columns
        ## these columns should have class character, the rest should have class numeric.
        # it only matters if I'm going to use bind_rows
        character_cols <- c("genome", "query", "hit", "desc_h", "strand_forGR")
        for (x in colnames(dat)) {
            if(x %in% character_cols) {
                dat[,x] <- as.character(dat[,x])
            } else {
                dat[,x] <- as.numeric(dat[,x])
            }
        }
        return(dat)
    }
    
    # convert strand to +/- so it's compatible with GRanges
    dat$strand_forGR <- NA
    dat[which(dat$strnd_h==1),"strand_forGR"] <- "+"
    dat[which(dat$strnd_h==-1),"strand_forGR"] <- "-"
    return(dat)
}


##### parsed_blast_output_to_gr = a function to make a GRanges for blast hits, parsed by the readBlast function. 
# Usually we want to get GRanges on the SUBJECT, not the QUERY, but I might provide the option later
parsed_blast_output_to_gr <- function(parsed_blast, 
                                      focus="subject") {
    if (focus=="subject") {
        
        # get seqlengths
        seqlen <- parsed_blast |> 
            select(seqnames=hit,
                   seqlengths=len_h) |> 
            unique() |> 
            as.data.frame() 
        seqlen <- Seqinfo(seqnames=seqlen$seqnames,
                          seqlengths=seqlen$seqlengths)
        
        # get hit regions
        parsed_blast_gr <- parsed_blast |> 
            dplyr::rename(seqnames=hit,
                          start=start_h, 
                          end=end_h,
                          strand=strand_forGR, 
                          identical=fr_idq,
                          conserved=fr_cnq) |> 
            select(-strnd_h, -strnd_q) |> 
            GRanges(seqinfo=seqlen)
        return(parsed_blast_gr)
    }
    if (focus=="query") {
        stop("\n\nERROR - haven't yet written code for this option\n\n")
    }
}

