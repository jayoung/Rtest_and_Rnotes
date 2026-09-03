library(ggmsa)

### xxx can I try this with devel version?

#### example protein alignment
protein_sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")

positions_to_highlight <- c(310,314,318)

#### ggmsa with position_highlight - this looks fine
ggmsa(protein_sequences, 
      start=300, end=350, 
      position_highlight = positions_to_highlight  ) 

#### add geom_msaBar - this looks fine
ggmsa(protein_sequences, 
      start=300, end=350) +
    geom_msaBar()

#### When we add position_highlight, two things happen that seems undesirable
# 1. the alignment view zooms in to the first/last highlighted position
# 2. the only columns shown in the geom_msaBar are the highlighted positions. Better to show all positions
ggmsa(protein_sequences, 
      start=300, end=350, 
      position_highlight = positions_to_highlight ) +
    geom_msaBar()

#### show sessionInfo()
sessionInfo()


