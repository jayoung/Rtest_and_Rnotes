## Issue 2: inconsistent scores on the diagonal

# https://github.com/Bioconductor/pwalign/issues/15

library(Biostrings)
library(pwalign)
data(BLOSUM62)


### Make small test alignment
aln_no_gaps <-  AAStringSet(c(seq1="VF",
                              seq2="VF",
                              seq3="LF"))

# seq1 and seq2 are identical to each other. 
# Therefore we should get identical scores for (a) seq1-versus-seq1 and (b) seq1-versus-seq2 .
# But those scores are not identical - self-matches (the diagonal) always get a score of 0.
# This makes sense for many distance metrics but maybe not for substitutionMatrix methods.

stringDist(aln_no_gaps, 
           diag=TRUE, upper=TRUE,
           method="substitutionMatrix",
           substitutionMatrix=BLOSUM62)

### Just to make sure I understand how the off-diagonal scores come about, here are the relevant bits of the BLOSUM62 matrix:
BLOSUM62[c("F","L","V"),c("F","L","V")]

### Understanding the scores:
# - VF aligned to VF scores -10, because BLOSUM62 scores F-F as 6 and V-V as 4 
# - VF aligned to LF score -7, because BLOSUM62 scores F-F as 6 and V-L as 1

### show sessionInfo
sessionInfo()

