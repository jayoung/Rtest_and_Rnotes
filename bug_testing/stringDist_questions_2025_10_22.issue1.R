## Issue 1: can StringDist handle gaps?

# https://github.com/Bioconductor/pwalign/issues/14

## Setup 
library(Biostrings)
library(pwalign)
data(BLOSUM62)

## Make small test alignments. `x` has no gaps, `y` has a gap.
aln_no_gaps <-  AAStringSet(c(seq1="VF",
                              seq2="VF",
                              seq3="LF"))

aln_with_gap <-  AAStringSet(c(seq1="VF",
                               seq2="V-",
                               seq3="LF"))

# stringDist works fine on the gapped alignment with default method (levenshtein)
stringDist(aln_with_gap,
           diag=TRUE, upper=TRUE)


# substitutionMatrix method does work without gaps:
stringDist(aln_no_gaps, 
           diag=TRUE, upper=TRUE,
           method="substitutionMatrix",
           substitutionMatrix=BLOSUM62)



# But substitutionMatrix method gives an error when there's a gap. Maybe the gap penalties are not being handled:
stringDist(aln_with_gap, 
           diag=TRUE, upper=TRUE,
           method="substitutionMatrix",
           substitutionMatrix=BLOSUM62)


# Finished
sessionInfo()

