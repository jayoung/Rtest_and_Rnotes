filename_processing_demo
================
Janet Young

2024-10-09

# Goal

Expore code given
[here](https://luisdva.github.io/rstats/filenames-to-vars/) that deals
with information contained within filenames (rather than within the
files)

Most of what’s below is just copied from that site.

# Group the data

Now, we can group the data by country of origin, year, and drive train
and then a) split the data into a list of tibbles (one for each group),
and b) use the group keys to build a vector of file names based on the
grouping information used to split the data.

``` r
# list of tibbles, 17 things
gtcars_groups <- gtcars |> 
    group_split(ctry_origin, year, drivetrain)

# tibble of the grouping variables
gtgroups <- 
    gtcars |> 
    group_by(ctry_origin, year, drivetrain) |> 
    # group_keys() returns a data frame describing the groups.
    group_keys() 
```

After minor changes like replacing spaces, the vector of filenames can
be built with rowwise glueing of the values in the grouping variables.

``` r
grpnames <- 
    gtgroups |> 
    rowwise() |> 
    mutate(gluedvars = glue_collapse(across(everything(), as.character), 
                                     sep = "_")) |> 
    mutate(gluedvars= str_replace(gluedvars," ","-")) |> 
    pull(gluedvars)
```

Let’s use this vector of filenames to name the list of tibbles, and then
export each tibble to a separate csv file using walk2 for its side
effects.

``` r
names(gtcars_groups) <- grpnames
gtcars_groups <- 
    gtcars_groups |> 
    map(select,-c(ctry_origin,year,drivetrain))
```

``` r
## walk2 is a purrr function - iterates over two variables at a time.
# save each to disk  
if(!dir.exists("filename_processing_demo_files")) {
    dir.create("filename_processing_demo_files")
}

walk2(gtcars_groups, # first thing to iterate over
      names(gtcars_groups), # second thing to iterate over
      # I assume the \ backslash is shorthand for function
      \(x,y) write_csv(x,paste0("filename_processing_demo_files/",y,".csv")))

## also create a vector of filenames
allpaths <- list.files("filename_processing_demo_files", 
                       full.names = TRUE)
```

Now let’s write a function that will strip the filename, split it into
fragments (one for each variable), and add these fragments as values
into the data rectangle as part of the import step. Using read.csv here
to get around guessing and type conversion.

``` r
filenameTovars <- function(filepath, var_names = NULL) {
    # strip the filename
    filename <- basename(filepath)
    filename <- tools::file_path_sans_ext(filename)
    
    # split the filename 
    filename_parts <- str_split(filename, "_", simplify = TRUE)
    
    # import
    data <- read_csv(filepath,
                     show_col_types = FALSE)
    
    # if no var_names provided
    if (is.null(var_names))  {
        var_names <- paste0("v", seq_along(filename_parts))
    }
    
    # filename fragments as new variables
    as_tibble(bind_cols(
        data,
        ## setNames() and unname() might be useful in piped context
        setNames(as.list(filename_parts), var_names) # a df that will be recycled
    ))
}
```

Nice little function above, it can take an optional vector with the
names of the variables being added and if none are provided, it will use
“v1”, “v2”, etc. Hard-coded “\_” as separator here, but that could be an
argument too.

For one of the files chosen at random, let’s compare the output from
reading the csv normally vs applying the new function with and without
the vector of variable names.

``` r
# read the csv normally
# read_csv("filename_processing_demo_files/Italy_2015_rwd.csv")

# read with path only (this puts the filename-derived variables in three columns called v1 v2 v3)
x <- filenameTovars("filename_processing_demo_files/Italy_2015_rwd.csv") 

# read with names for new columns
x <- filenameTovars("filename_processing_demo_files/Italy_2015_rwd.csv",
                    var_names = c("ctry_origin","year","drivetrain")) 
```

Finally, we can use purrr to iterate through all the files in a folder
and apply the custom function. Note the new recommended approach of
using list_rbind rather than map_df.

``` r
dat <- map(allpaths,
           \(x) 
           filenameTovars(x,
                          var_names = c("ctry_origin","year",
                                        "drivetrain")))  |> 
    ## I needed to add this as.character() line, otherwise the list_rbind failed
    # Can't combine `..1$trim` <character> and `..16$trim` <double>.
    map(\(x) mutate(x, trim=as.character(trim)) ) |>
    list_rbind() |>
    ## type_convert retries the process of guessing what class each column of the data should be
    suppressMessages(type_convert())
```

# Tidy up

``` r
dir_delete("filename_processing_demo_files")
```

# Finished

show R version used, and package versions

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: x86_64-apple-darwin20
    ## Running under: macOS Ventura 13.7
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] fs_1.6.4      readr_2.1.5   stringr_1.5.1 glue_1.8.0    purrr_1.0.2  
    ## [6] dplyr_1.1.4   gt_0.11.1    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] crayon_1.5.3      vctrs_0.6.5       cli_3.6.3         knitr_1.47       
    ##  [5] rlang_1.1.4       xfun_0.45         stringi_1.8.4     generics_0.1.3   
    ##  [9] bit_4.0.5         htmltools_0.5.8.1 hms_1.1.3         fansi_1.0.6      
    ## [13] rmarkdown_2.27    evaluate_0.24.0   tibble_3.2.1      tzdb_0.4.0       
    ## [17] fastmap_1.2.0     yaml_2.3.8        lifecycle_1.0.4   compiler_4.4.0   
    ## [21] pkgconfig_2.0.3   rstudioapi_0.16.0 digest_0.6.36     R6_2.5.1         
    ## [25] tidyselect_1.2.1  utf8_1.2.4        parallel_4.4.0    vroom_1.6.5      
    ## [29] pillar_1.9.0      magrittr_2.0.3    bit64_4.0.5       withr_3.0.1      
    ## [33] tools_4.4.0       xml2_1.3.6
