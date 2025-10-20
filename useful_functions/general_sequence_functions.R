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
        seqs <- seq_file %>% readDNAStringSet(nrec=nrec)
    }
    if (type=="AA") {
        seqs <- seq_file %>% readAAStringSet(nrec=nrec)
    }
    if (type=="B") {
        seqs <- seq_file %>% readBStringSet(nrec=nrec)
    }
    if(!exists("seqs")) {
        stop("\n\nERROR - you specified an invalid type. Should be DNA, AA, or B\n\n")
    }
    headers_split <- names(seqs) %>% strsplit(" ")
    info_tbl <- tibble(seq_id = sapply(headers_split, "[[", 1),
                       seq_len = width(seqs),
                       desc = sapply(headers_split, function(x) {
                           if(length(x)>1) {
                               return( paste(x[2:length(x)], collapse=" ")  )
                           } else { return(NA) }
                       }))
    if(strip_afterDot) {
        info_tbl <- info_tbl %>% 
            mutate(seq_id_withDot = seq_id) %>% 
            mutate(seq_id = str_remove_all(seq_id, "\\.\\d+$")) %>% 
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
# test[["seqs"]] %>% names()
# test[["info"]]


#### countAmbiguitiesEachSeq - a small function that counts all non-ACGT seqs in each member of a DNAStringSet
countAmbiguitiesEachSeq <- function(dna_seqs) {
    alphabetFrequency(dna_seqs) %>% 
        as_tibble() %>% 
        select(-A, -C, -G, -T) %>% 
        rowSums()
}


