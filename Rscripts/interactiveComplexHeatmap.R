# https://github.com/jokergoo/InteractiveComplexHeatmap
# https://jokergoo.github.io/InteractiveComplexHeatmap/articles/deseq2_app.html

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("airway")

library(InteractiveComplexHeatmap)
library(DESeq2)
library(airway)

data(airway)
dds = DESeqDataSet(airway, design = ~ dex)
dds = DESeq(dds)

interactivate(dds)




### design: 
# colData(airway)
# 8 samples, there's a 'dex' column - 4 'untrt' and 4 'trt' samples

### counts has 64102 rows and 8 columns
# dim(assay(airway))
