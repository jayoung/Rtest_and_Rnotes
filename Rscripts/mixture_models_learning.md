mixture_models_learning
================
Janet Young

2026-05-19

# Goal

Learn about mixture modelling and explore related R packages.

<https://en.wikipedia.org/wiki/Mixture_model>

<https://yifengedms.github.io/EDMS657-R-Tutorials/Mixture.html>

<https://mclust-org.github.io/mclust-book/>

# General notes

Often we model a mixture of subpopulations that all follow the same type
of distribution (e.g. all normal) but it is also possible to model
mixtures of different types of distribution.

Expectation maximization algorithms are often used to do the modeling

See also [a
script](https://github.com/gak2882/SATAY/blob/main/260615_flowData_analysis/Rscripts_janet/100_flowData_analysis_janet.Rmd)
and [its
output](https://github.com/gak2882/SATAY/blob/main/260615_flowData_analysis/Rscripts_janet/100_flowData_analysis_janet.Rmd)
that I wrote to looks at flow sorting data for the 2micron SATAY project

# Packages used

mclust
[vignette](https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html)

mixtools
[vignette](https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf)
and [github](https://github.com/dsy109/mixtools) page

[plotmm
package](https://packages.oit.ncsu.edu/cran/web/packages/plotmm/vignettes/Getting-Started.html)
for plotting model results

# Example data

The classic example dataset is the wait times between eruptions of the
Old Faithful geyser, which seems to follow a bimodal distribution. See
`?faithful` for more information.

Show the data we are modelling:

``` r
wait_histo <- faithful |> 
    ggplot(aes(x=waiting)) +
    geom_histogram(breaks=seq(from=40, to=100, by=5),
                   fill="lightgray", color="black", linewidth=0.2) +
    theme_classic() +
    labs(x="Wait time (minutes)",
         y="number of observations",
         title="(histogram)")

wait_density <- faithful |> 
    ggplot(aes(x=waiting)) +
    geom_density() +
    theme_classic() +
    labs(x="Wait time (minutes)",
         y="Relative frequency",
         title="density plot")

(wait_histo + wait_density) +
    plot_annotation(title="Old Faithful data,\ndistribution of wait time between eruptions")
```

![](mixture_models_learning_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

# Try mclust package

Use the `densityMclust` function to model the data. The default setup
(for univariate data, like we have) is:

- it returns a `densityMclust` class object
- it models different numbers of classes (components), from 1-9
- it models using equal “E” or unequal (“V”) variance
- it can produces a bunch of other plot types, see `?plot.densityMclust`

``` r
faithful_mclust <- densityMclust(faithful$waiting,
                                 plot=FALSE, # default is to plot the density of the data
                                 verbose=FALSE)
```

`summary()` is inelegant but shows me the best model:

``` r
summary(faithful_mclust, parameters = TRUE)
```

    ## ------------------------------------------------------- 
    ## Density estimation via Gaussian finite mixture modeling 
    ## ------------------------------------------------------- 
    ## 
    ## Mclust E (univariate, equal variance) model with 2 components: 
    ## 
    ##  log-likelihood   n df       BIC       ICL
    ##       -1034.002 272  4 -2090.427 -2099.576
    ## 
    ## Mixing probabilities:
    ##         1         2 
    ## 0.3609461 0.6390539 
    ## 
    ## Means:
    ##        1        2 
    ## 54.61675 80.09239 
    ## 
    ## Variances:
    ##        1        2 
    ## 34.44093 34.44093

We can also ask it to plot the BIC of each model it considered.

Here we see that two classes is the the most likely solution fit, and
that equal variance is slightly more likely than unequal variance.

BIC is a metric that includes some penalty for each additional parameter
in the model, so that it tries to avoid overfitting.

We can show the model and the underlying data together:

``` r
par(mfrow=c(1,2))
plot(faithful_mclust, what = "BIC")
title(main="model selection", outer=FALSE, line=0.5)

plot(faithful_mclust, what = "density", data = faithful$waiting)
title(main="best model", outer=FALSE, line=0.5)

title(main="mclust modelling", outer=TRUE, line=-2)
```

![](mixture_models_learning_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

“Diagnostic” plots show actual and modeled distributions.

``` r
par(mfrow=c(1,2))

plot(faithful_mclust, what = "diagnostic", type="cdf")
title(main="CDF diagnostic plot", line=0.5)

plot(faithful_mclust, what = "diagnostic", type="qq")
title(main="QQ diagnostic plot", line=0.5)
```

![](mixture_models_learning_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Made a few attempts to plot the individual components of an mclust model
but I don’t have it working yet. I bet there’s some package that does it
that I haven’t found yet.

``` r
## verbose FALSE so we don't get the progress message
faithful_mclust_G2 <- densityMclust(faithful$waiting, G = 2, plot = FALSE, verbose=FALSE)
```

``` r
### here's how I can extract the mean and standard deviation of each component and plot them
faithful_x_range_for_density <- extendrange(faithful$waiting, f = 0.3)
faithful_x_values_for_density <- seq(from=faithful_x_range_for_density[1], 
                                     to=faithful_x_range_for_density[2], 
                                     length=1000)

y1 <- dnorm(faithful_x_values_for_density, 
            mean=faithful_mclust_G2$parameters$mean[1], 
            sd=sqrt(faithful_mclust_G2$parameters$variance$sigmasq)) * faithful_mclust_G2$parameters$pro[1]  

y2 <- dnorm(faithful_x_values_for_density, 
            mean=faithful_mclust_G2$parameters$mean[2], 
            sd=sqrt(faithful_mclust_G2$parameters$variance$sigmasq)) * faithful_mclust_G2$parameters$pro[2]  

faithful_mclust_G2_model_to_plot <- tibble(x=faithful_x_values_for_density, 
       y1=y1, 
       y2=y2) |> 
    mutate(y_total=y1+y2)

faithful_mclust_G2_model_to_plot |> 
    ggplot(aes(x=x, y=y1)) +
    geom_histogram(data=faithful, aes(x=waiting, y=after_stat(density)), 
                   fill="lightgray", color="darkgray",
                   binwidth = 5) +
    geom_density(data=faithful, aes(x=waiting, y=after_stat(density)), 
                   color="darkgray") +
    geom_line(color="red") +
    geom_line(aes(x=x, y=y2), color="blue") +
    geom_line(aes(x=x, y=y_total), color="purple") +
    theme_classic() +
    labs(x="Waiting time (minutes)",
         y="Relative frequency",
         title="Old Faithful waiting data and mclust G2 model",
         subtitle=paste0("Gray histogram and density = real data\n",
                         "Red/blue/purple = modeled components"))
```

![](mixture_models_learning_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# mclust1Dplot is only for multivariate models
# mclust1Dplot(faithful_mclust_G2)
```

# Try mixtools package

Do the modelling using `mixtools::normalmixEM`:

- the result (`wait1`) is a `mixEM` object
- `mu` - we provide initial estimates for the means of each
  distribution. Because we start with a vector of 2 for mu, it models a
  mix of 2 normal distributions
- `lambda` is the initial mixing proportion. (if you don’t specify,
  it’ll start with equal shares)
- it does 9 interations as it optimizes the proportions and means and
  sigma
- `sigma` is the starting standard deviation
- you can use `mean.constr` to constrain one or more of the means, which
  could be useful in our case, where we could use the cir0 distribution
  to guess at one of the components. Same for `sd.constr`

Let’s tell it there are 2 categories, but nothing else (normalmixEM does
not allow us to try different k all in one call, like mclust does). By
default, the variance to be different on the two components but you can
constrain it to be the same. You can specify starting mu (mean) and
sigma (standard deviation), and you can constrain some components but
not others (using mean.constr and sd.constr, see ?normalmixEM)

Show a summary of the model:

``` r
summary(wait1)
```

    ## summary of normalmixEM object:
    ##          comp 1   comp 2
    ## lambda  0.36085  0.63915
    ## mu     54.61364 80.09031
    ## sigma   5.86909  5.86909
    ## loglik at estimate:  -1034.002

Plot the model and the input data:

- plot() calls plot.mixEM (`?plot.mixEM`)
- by default this makes two plots (density and loglik) and requires the
  user to hit return between each, but we can control that:
- density plot (left) shows the component distributions
- `loglik` plot (right) shows how the log-likelihood of the model
  changed as the ME iterated - after 2 iterations it didn’t improve much

``` r
par(mfrow=c(1,2))
plot(wait1, 
     loglik=FALSE, density=TRUE,
     main2="Model (curves) and data (histogram)" ,
     xlab2="Minutes")

plot(wait1, 
     loglik=TRUE, density=FALSE, main1="Log-likelihoods")
```

![](mixture_models_learning_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

plotmm::plot_mm can also plot the model (uses ggplot approach). the gray
shape is the distribution of the data being modelled, and the colors are
the actual data

Works on mixEM objects (e.g. output of mixtools::normalmixEM). Also
supports objects modeled by ‘EMCluster’, and ‘flexmix’ packages (but not
mclust).

``` r
p1 <- plot_mm(wait1, 2) +
    labs(title = "normalmixEM model for faithful$waiting",
         subtitle = "Plotted using plot_mm") +
    coord_cartesian(xlim=c(40,100)) 

##  I checked - the gray shows density of the data that was modeled:
# p1 +
#   geom_density(data=faithful, aes(x=waiting), 
#                color="forestgreen", lty=2)

p1
```

![](mixture_models_learning_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Can also customize that using `plot_mix_comps_normal()` as follows:

``` r
data.frame(x = wait2$x) |>
    ggplot() +
    ## show distribution of the actual data
    geom_histogram(aes(x, after_stat(density)), 
                   binwidth = 5, colour = "black", fill="lightgray") +
    ### plot component 1
    stat_function(geom = "line", 
                  fun = plot_mix_comps_normal, # here is the function
                  args = list(wait2$mu[1], wait2$sigma[1], lam = wait2$lambda[1]),
                  colour = "red", lwd = 1.5) +
    ### plot component 2
    stat_function(geom = "line", 
                  fun = plot_mix_comps_normal, # here again as k = 2
                  args = list(wait2$mu[2], wait2$sigma[2], lam = wait2$lambda[2]),
                  colour = "blue", lwd = 1.5) +
    ylab("Density") +
    theme_classic()
```

![](mixture_models_learning_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

`plot_cut_point()` - if we want to use the model to predict which class
each datapoint is in, we might want to determine the thresholds between
classes:

``` r
## message=FALSE otherwise I get a message about the binwidth
plot_cut_point(wait1, plot = TRUE, color = "amerika") + # produces plot
    labs(x="Wait time (minutes)")
```

![](mixture_models_learning_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

mixtools has a function called multmixmodel.sel that assesses how many
components are in the data, but that’s only for multivariate data. I
don’t see an equivalent function for univariate data, unless we have \>1
sample of the data. I think maybe we are supposed to subsample for a
bootstrapping approach, see ?boot.se

### mixtools::mixturegram

See ?mixturegram

We generate some example data that’s a 30:70 mix of two normal
distributions, with means of -6 and 0, and standard deviation of 1 for
both

``` r
set.seed(100)
n <- 100
w <- rmultinom(n,1,c(.3,.7))
y <- sapply(1:n, function(i)  {
    w[1,i]*rnorm(1,-6,1) + w[2,i]*rnorm(1,0,1)
})
```

Show the distribution of the data:

``` r
tibble(y=y) |> 
    ggplot(aes(x=y)) +
    geom_density() +
    theme_classic()
```

![](mixture_models_learning_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Show a ‘mixturegram’. I don’t understand what the “PC score” is here.

``` r
selection <- function(i, data, rep=30){
    out <- replicate(rep,normalmixEM(data,epsilon=1e-06,
                                     k=i,maxit=5000),simplify=FALSE)
    counts <- lapply(1:rep,function(j) 
        table(apply(out[[j]]$posterior,1,
                    which.max)))
    counts.length <- sapply(counts, length)
    counts.min <- sapply(counts, min)
    counts.test <- (counts.length != i)|(counts.min < 5)
    if(sum(counts.test) > 0 & sum(counts.test) < rep) 
        out <- out[!counts.test]
    l <- unlist(lapply(out, function(x) x$loglik))
    tmp <- out[[which.max(l)]]
}

### does mixture modelling with 2:5 classes, replicating each 30 times
all.out <- lapply(2:5, selection, data = y, rep = 2)
```

    ## number of iterations= 7 
    ## number of iterations= 7 
    ## number of iterations= 52 
    ## number of iterations= 88 
    ## number of iterations= 266 
    ## number of iterations= 111 
    ## number of iterations= 262 
    ## number of iterations= 69

``` r
## pmbs is a list of length 4, each is a matrix, 100 rows by 2,3,4,5 columns
pmbs <- lapply(1:length(all.out), function(i) 
    all.out[[i]]$post)


mixturegram(y, pmbs, method = "pca", all.n = FALSE,
            id.con = NULL, score = 1, 
            main = "Mixturegram (Well-Separated Data)")
```

![](mixture_models_learning_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

    ## $stopping
    ## [1] 1.000000000 0.012826548 0.005056038 0.004510841 0.003132543

# Try tidyLPA

I’m trying the `estimate_profiles` function from the `tidyLPA` package,
which is a wrapper for mclust’s mixture models

See
[tutorial](https://yifengedms.github.io/EDMS657-R-Tutorials/Mixture.html#Example_1:_Univariate_Mixtures)
for a very brief outline.

See [tidyLPA
intro](https://cran.r-project.org/web/packages/tidyLPA/vignettes/Introduction_to_tidyLPA.html)

``` r
faithful_tidyLPA <- estimate_profiles(faithful$waiting, 
                              ## this says, try it with 1-6 mixture components
                              n_profiles=1:6,
                              variances = c("equal", "varying"),
                              covariances = c("zero", "zero"))
## by default it says it does the modelling using 'OpenMx' package but it can also use mclust or MplusAutomation. ?estimate_profiles is a bit ambiguous about what the default is, actually.
```

Here’s the overall result:

``` r
faithful_tidyLPA
```

    ## tidyLPA analysis using mclust: 
    ## 
    ##    Model Classes     AIC     BIC Entropy prob_min prob_max n_min n_max BLRT_p
    ## 1      1       1 2194.58 2201.79    1.00     1.00     1.00  1.00  1.00       
    ## 2      1       2 2076.00 2090.43    0.94     0.98     0.99  0.36  0.64   0.01
    ## 3      1       3 2080.21 2101.85    0.59     0.43     0.99  0.21  0.43   0.86
    ## 4      1       4 2079.50 2108.35    0.72     0.53     0.93  0.11  0.50   0.16
    ## 5      1       5 2082.40 2118.46    0.67     0.50     0.93  0.08  0.42   0.24
    ## 6      1       6 2086.50 2129.77    0.61     0.35     0.93  0.05  0.33   0.78
    ## 7      2       1 2194.58 2201.79    1.00     1.00     1.00  1.00  1.00       
    ## 8      2       2 2078.01 2096.04    0.94     0.98     0.99  0.36  0.64   0.01
    ## 9      2       3 2084.15 2113.00    0.66     0.70     0.98  0.31  0.37   0.93
    ## 10     2       4 2087.05 2126.72    0.57     0.53     0.92  0.13  0.33   0.19
    ## 11     2       5 2090.33 2140.82    0.64     0.64     0.91  0.15  0.26   0.31
    ## 12     2       6 2091.53 2152.83    0.69     0.60     0.89  0.12  0.20   0.24

``` r
faithful_tidyLPA |> 
    compare_solutions()
```

    ## Compare tidyLPA solutions:
    ## 
    ##  Model Classes BIC     
    ##  1     1       2201.789
    ##  1     2       2090.427
    ##  1     3       2101.849
    ##  1     4       2108.345
    ##  1     5       2118.459
    ##  1     6       2129.772
    ##  2     1       2201.789
    ##  2     2       2096.044
    ##  2     3       2112.995
    ##  2     4       2126.717
    ##  2     5       2140.816
    ##  2     6       2152.833
    ## 
    ## Best model according to BIC is Model 1 with 2 classes.
    ## 
    ## An analytic hierarchy process, based on the fit indices AIC, AWE, BIC, CLC, and KIC (Akogul & Erisoglu, 2017), suggests the best solution is Model 1 with 2 classes.

`get_data()` gets the class assignments for each datapoint for each
model. Here I choose the best model, and the output has 544 rows (two
per observation in faithful\$waiting - probability of being in each
class)

``` r
faithful_tidyLPA |> 
    get_data() |> 
    filter(model_number==1 & classes_number==2) |> 
    head()
```

    ## # A tibble: 6 × 7
    ##   model_number classes_number    df Class Class_prob Probability    id
    ##          <dbl>          <int> <dbl> <dbl>      <dbl>       <dbl> <int>
    ## 1            1              2    79     2          1  0.000103       1
    ## 2            1              2    54     1          1  1.000          2
    ## 3            1              2    74     2          1  0.00415        3
    ## 4            1              2    62     1          1  0.968          4
    ## 5            1              2    85     2          1  0.00000122     5
    ## 6            1              2    55     1          1  1.000          6

``` r
p1 <- faithful_tidyLPA |> 
    get_fit() |> 
    mutate(Model=as.factor(Model)) |> 
    ggplot(aes(x=Classes, y=(0-BIC), #group=Model, 
               color=Model)) +
    geom_line() +
    geom_point() +
    theme_classic()+
    labs(title="BIC")
p2 <- faithful_tidyLPA |> 
    get_fit() |> 
    mutate(Model=as.factor(Model)) |> 
    ggplot(aes(x=Classes, y=(0-AIC), #group=Model, 
               color=Model)) +
    geom_line() +
    geom_point() +
    theme_classic() +
    labs(title="AIC")
(p1 + p2) +
    plot_layout(guides="collect")
```

![](mixture_models_learning_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
faithful_tidyLPA |> 
    get_estimates() |> 
    filter(Model==1 & Classes==2) |> 
    select(-Parameter, -Model, -Classes) |> 
    relocate(Class)
```

    ## # A tibble: 4 × 5
    ##   Class Category  Estimate    se        p
    ##   <int> <chr>        <dbl> <dbl>    <dbl>
    ## 1     1 Means         54.6 0.602 0       
    ## 2     1 Variances     34.4 3.18  2.68e-27
    ## 3     2 Means         80.1 0.526 0       
    ## 4     2 Variances     34.4 3.18  2.68e-27

``` r
# faithful_mclust_G2 |> 
     # filter(Model==1 & Classes==2) |> 
    # plot_profiles()
    # plot_density()
```

# TidyDensity package

`tidy_normal()` might be helpful.

See [blog
post](https://www.spsanderson.com/steveondata/posts/2025-06-23/) or
[vignette](https://cran.r-project.org/web/packages/TidyDensity/vignettes/getting-started.html)

``` r
tidy_normal(.n = 200, 
            .mean = 0, .sd = 1, 
            .num_sims = 1, .return_tibble = TRUE) |> 
    ggplot(aes(x=y)) +
    geom_density() +
    theme_classic()
```

![](mixture_models_learning_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

Use tidy_normal to get a mixture distribution:

``` r
dat <- bind_rows(
    tidy_normal(.n = 100, 
                   .mean = 0, .sd = 1, 
                   .num_sims = 1, .return_tibble = TRUE) |> 
        mutate(my_sim_num="A"),
    tidy_normal(.n = 500, 
                   .mean = 10, .sd = 1, 
                   .num_sims = 1, .return_tibble = TRUE)  |> 
        mutate(my_sim_num="B"))

## geom_density default is to scaling so that each group has the same area, 
# use after_stat(count) to scale accordingly
dat |>
    ggplot(aes(x=y, y=after_stat(count), group=my_sim_num, color=my_sim_num)) +
    geom_density(position="stack") +
    theme_classic()
```

![](mixture_models_learning_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

# Finished

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Tahoe 26.5
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
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
    ##  [1] TidyDensity_1.5.2 tidyLPA_2.0.2     tidySEM_0.2.10    mixtools_2.0.0.1 
    ##  [5] plotmm_0.1.2      mclust_6.1.2      patchwork_1.3.2   here_1.0.2       
    ##  [9] kableExtra_1.4.0  lubridate_1.9.5   forcats_1.0.1     stringr_1.6.0    
    ## [13] dplyr_1.2.1       purrr_1.2.1       readr_2.2.0       tidyr_1.3.2      
    ## [17] tibble_3.3.1      ggplot2_4.0.2     tidyverse_2.0.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] mnormt_2.1.2        sandwich_3.1-1      rlang_1.1.7        
    ##   [4] magrittr_2.0.5      otel_0.2.0          compiler_4.5.3     
    ##   [7] flexmix_2.3-20      systemfonts_1.3.2   vctrs_0.7.2        
    ##  [10] quadprog_1.5-8      crayon_1.5.3        pkgconfig_2.0.3    
    ##  [13] fastmap_1.2.0       backports_1.5.1     labeling_0.4.3     
    ##  [16] utf8_1.2.6          pbivnorm_0.6.0      pander_0.6.6       
    ##  [19] rmarkdown_2.31      tzdb_0.5.0          modeltools_0.2-24  
    ##  [22] xfun_0.57           jsonlite_2.0.0      progress_1.2.3     
    ##  [25] EMCluster_0.2-17    psych_2.6.5         prettyunits_1.2.0  
    ##  [28] parallel_4.5.3      lavaan_0.6-21       R6_2.6.1           
    ##  [31] stringi_1.8.7       RColorBrewer_1.1-3  parallelly_1.46.1  
    ##  [34] car_3.1-5           boot_1.3-32         Rcpp_1.1.1         
    ##  [37] knitr_1.51          future.apply_1.20.2 zoo_1.8-15         
    ##  [40] nnet_7.3-20         Matrix_1.7-5        splines_4.5.3      
    ##  [43] igraph_2.2.2        timechange_0.4.0    tidyselect_1.2.1   
    ##  [46] rstudioapi_0.18.0   abind_1.4-8         yaml_2.3.12        
    ##  [49] codetools_0.2-20    listenv_0.10.1      lattice_0.22-9     
    ##  [52] nonnest2_0.5-9      plyr_1.8.9          withr_3.0.2        
    ##  [55] S7_0.2.1            coda_0.19-4.1       evaluate_1.0.5     
    ##  [58] future_1.70.0       fastDummies_1.7.6   survival_3.8-6     
    ##  [61] CompQuadForm_1.4.4  xml2_1.5.2          texreg_1.39.5      
    ##  [64] kernlab_0.9-33      pillar_1.11.1       carData_3.0-6      
    ##  [67] checkmate_2.3.4     stats4_4.5.3        plotly_4.12.0      
    ##  [70] generics_0.1.4      dbscan_1.2.4        rprojroot_2.1.1    
    ##  [73] hms_1.1.4           scales_1.4.0        globals_0.19.1     
    ##  [76] xtable_1.8-8        glue_1.8.0          lazyeval_0.2.3     
    ##  [79] tools_4.5.3         data.table_1.18.2.1 amerika_0.1.1      
    ##  [82] gsubfn_0.7          RANN_2.6.2          mvtnorm_1.3-7      
    ##  [85] grid_4.5.3          MplusAutomation_1.2 wesanderson_0.3.7  
    ##  [88] nlme_3.1-169        proto_1.0.0         Formula_1.2-5      
    ##  [91] cli_3.6.5           textshaping_1.0.5   segmented_2.2-1    
    ##  [94] viridisLite_0.4.3   svglite_2.2.2       gtable_0.3.6       
    ##  [97] digest_0.6.39       progressr_0.19.0    htmlwidgets_1.6.4  
    ## [100] farver_2.1.2        htmltools_0.5.9     lifecycle_1.0.5    
    ## [103] httr_1.4.8          MASS_7.3-65
