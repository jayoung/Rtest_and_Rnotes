# Load packages 
library("tidyverse")
library("patchwork")

# Make some plots
dat <- tibble(x=1:100, y=1:100)

p1 <- dat |> 
    ggplot(aes(x=x, y=y)) +
    geom_point() + 
    labs(title="p1")

p2 <- p1 + labs(title="p2")
p3 <- p1 + labs(title="p3")

# Show plots in desired layout - works fine with only tidyverse and patchwork loaded
(p1 / p2) | p3  

# Use `+` not `|` and we get a different layout (no error)
(p1 / p2) + p3  

# Now we load ggedit package (but don't edit the plots at all)
library("ggedit") 

# Now patchwork fails when we use `|` rather than `+`
(p1 / p2) | p3 # , fig.height=5, fig.width=9
## Error: There are no null layers available in the plot to remove


# It works if we replace `|` with `+` this isn't the layout we want
(p1 / p2) + p3

# If we try using plot_layout (patchwork) we get a different error (even though we're still using | not +)
(p1 / p2) | p3 + plot_layout(widths=c(1,2))
# Error in `l[[x]]`:
#     ! Patchworks can only be indexed with numeric indices
# Run `rlang::last_trace()` to see where the error occurred.

# Show `sessionInfo()`
sessionInfo()
# R version 4.5.0 (2025-04-11)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.5
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
#     [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Los_Angeles
# tzcode source: internal
# 
# attached base packages:
#     [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#     [1] ggedit_0.4.1    patchwork_1.3.0 lubridate_1.9.4 forcats_1.0.0  
# [5] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.4     readr_2.1.5    
# [9] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.2   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
#     [1] jsonlite_2.0.0     miniUI_0.1.2       gtable_0.3.6      
# [4] compiler_4.5.0     promises_1.3.3     Rcpp_1.0.14       
# [7] tidyselect_1.2.1   shinyAce_0.4.4     later_1.4.2       
# [10] scales_1.4.0       fastmap_1.2.0      mime_0.13         
# [13] R6_2.6.1           labeling_0.4.3     generics_0.1.4    
# [16] shiny_1.10.0       pillar_1.10.2      RColorBrewer_1.1-3
# [19] tzdb_0.5.0         rlang_1.1.6        stringi_1.8.7     
# [22] httpuv_1.6.16      timechange_0.3.0   cli_3.6.5         
# [25] withr_3.0.2        magrittr_2.0.3     shinyBS_0.61.1    
# [28] digest_0.6.37      grid_4.5.0         xtable_1.8-4      
# [31] rstudioapi_0.17.1  hms_1.1.3          lifecycle_1.0.4   
# [34] vctrs_0.6.5        glue_1.8.0         farver_2.1.2      
# [37] tools_4.5.0        pkgconfig_2.0.3    htmltools_0.5.8.1 
