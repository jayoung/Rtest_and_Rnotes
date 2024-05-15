Why use Rmarkdown?
================
Janet Young

2024-05-15

# What is Rmarkdown?

Rmarkdown (Rmd) helps us generate data analysis code that is
reproducible and well-documented.

With a Rmd document, we can run the code interactively (in the console).
Or, once our code is working, we can use “knit” to generate a report in
the format(s) of our choice - html, word doc, etc.

We can also easily choose to “restart R and run all chunks”. This
ensures the code commands are all in an order that makes sense, and that
the output will be reproducible in a “clean” environment. Running code
interactively runs the risk we’ve run code blocks in a different order,
or loaded packages or data that we haven’t specified in the script, etc.

We control the output report format using the YAML block at the top of
the script - I like to use the following output format, because when I
sync my project to github, the report displays nicely. I can also find
the files for individual plots within the report, in a subdirectory
nearby.

    output: github_document
    always_allow_html: true

# Alternatives to Rmarkdown

Some people like **Jupyter notebooks** - those can run R code (or
python)

**Quarto** is similar Rmarkdown - I’ve only used it briefly. Similar
concept to Rmd but formatting might be easier. It looks nice, but I
don’t know it it’s mature enough for regular use. It probably is. It’s
available in later versions of Rstudio.

# How to use Rmd

I’m not going to explain the whole thing here - any of these tutorials
should help you learn it fairly quickly:

- Rstudio’s [intro to
  Rmarkdown](https://rmarkdown.rstudio.com/lesson-1.html)  
- [intro2r chapter 8](https://intro2r.com/rmarkdown_r.html)  
- detailed [Rmarkdown guide](https://bookdown.org/yihui/rmarkdown/)

I just want to use this script to point out a few useful features.

We often start with a code chunk that sets “global” options and loads
packages we’ll need in the script. This `echo = TRUE` controls whether
your code itself will be visible in your output document. Sometimes it’s
nice to see it, sometimes you prefer not to.

``` r
knitr::opts_chunk$set(echo = TRUE)

## if you want to include code that doesn't need to be run every time, but maybe you want to remind yourself how you did something, comment it out, otherwise it will run every time you run the script, slowing things down. Common example is installing packages
# install.packages("tidyverse")

library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(kableExtra) # this package helps display tables nicely
```

    ## 
    ## Attaching package: 'kableExtra'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     group_rows

Giving each code chunk a short description (e.g. the `setup` in
`{r setup}` above) can help you see which code chunks are slow.

If the script is running too slow for my liking, I sometimes split up my
code into a couple of scripts - the first one runs any slow data
crunching, and saves the output. The second one loads the saved output
and makes plots/tables etc.

# Example

We’ll use R’s built-in dataset `iris` just to demo some things.

Stuff that would print onto your R console can print to your report too

``` r
head(iris)
```

    ##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    ## 1          5.1         3.5          1.4         0.2  setosa
    ## 2          4.9         3.0          1.4         0.2  setosa
    ## 3          4.7         3.2          1.3         0.2  setosa
    ## 4          4.6         3.1          1.5         0.2  setosa
    ## 5          5.0         3.6          1.4         0.2  setosa
    ## 6          5.4         3.9          1.7         0.4  setosa

``` r
class(iris)
```

    ## [1] "data.frame"

``` r
dim(iris)
```

    ## [1] 150   5

Displaying a table via kable is a bit prettier. There are also sorts of
ways to control formatting here, too - google will help with that

``` r
iris %>% 
    kable() %>% 
    kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
Sepal.Length
</th>
<th style="text-align:right;">
Sepal.Width
</th>
<th style="text-align:right;">
Petal.Length
</th>
<th style="text-align:right;">
Petal.Width
</th>
<th style="text-align:left;">
Species
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
3.5
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.9
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.7
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.6
</td>
<td style="text-align:right;">
3.1
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
3.6
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.4
</td>
<td style="text-align:right;">
3.9
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.6
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.4
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.9
</td>
<td style="text-align:right;">
3.1
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.4
</td>
<td style="text-align:right;">
3.7
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.8
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.8
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.3
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
1.1
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.8
</td>
<td style="text-align:right;">
4.0
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
4.4
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.4
</td>
<td style="text-align:right;">
3.9
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
3.5
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
3.8
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
3.8
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.4
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
3.7
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.6
</td>
<td style="text-align:right;">
3.6
</td>
<td style="text-align:right;">
1.0
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
3.3
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.8
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.2
</td>
<td style="text-align:right;">
3.5
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.2
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.7
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.8
</td>
<td style="text-align:right;">
3.1
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.4
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.2
</td>
<td style="text-align:right;">
4.1
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.5
</td>
<td style="text-align:right;">
4.2
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.9
</td>
<td style="text-align:right;">
3.1
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.5
</td>
<td style="text-align:right;">
3.5
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.9
</td>
<td style="text-align:right;">
3.6
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.4
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
3.5
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.4
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
3.5
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:right;">
0.6
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
3.8
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.8
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
3.8
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
4.6
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.3
</td>
<td style="text-align:right;">
3.7
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
3.3
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:left;">
setosa
</td>
</tr>
<tr>
<td style="text-align:right;">
7.0
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
4.7
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.4
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.9
</td>
<td style="text-align:right;">
3.1
</td>
<td style="text-align:right;">
4.9
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.5
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:right;">
4.0
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.5
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
4.6
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.3
</td>
<td style="text-align:right;">
3.3
</td>
<td style="text-align:right;">
4.7
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
4.9
</td>
<td style="text-align:right;">
2.4
</td>
<td style="text-align:right;">
3.3
</td>
<td style="text-align:right;">
1.0
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.6
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
4.6
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.2
</td>
<td style="text-align:right;">
2.7
</td>
<td style="text-align:right;">
3.9
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
2.0
</td>
<td style="text-align:right;">
3.5
</td>
<td style="text-align:right;">
1.0
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.9
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
4.2
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.0
</td>
<td style="text-align:right;">
2.2
</td>
<td style="text-align:right;">
4.0
</td>
<td style="text-align:right;">
1.0
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
4.7
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
3.6
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
3.1
</td>
<td style="text-align:right;">
4.4
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.8
</td>
<td style="text-align:right;">
2.7
</td>
<td style="text-align:right;">
4.1
</td>
<td style="text-align:right;">
1.0
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.2
</td>
<td style="text-align:right;">
2.2
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:right;">
3.9
</td>
<td style="text-align:right;">
1.1
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.9
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
4.8
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
4.0
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.3
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:right;">
4.9
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
4.7
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.4
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
4.3
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.6
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
4.4
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.8
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
4.8
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.0
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
2.6
</td>
<td style="text-align:right;">
3.5
</td>
<td style="text-align:right;">
1.0
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.5
</td>
<td style="text-align:right;">
2.4
</td>
<td style="text-align:right;">
3.8
</td>
<td style="text-align:right;">
1.1
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.5
</td>
<td style="text-align:right;">
2.4
</td>
<td style="text-align:right;">
3.7
</td>
<td style="text-align:right;">
1.0
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.8
</td>
<td style="text-align:right;">
2.7
</td>
<td style="text-align:right;">
3.9
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.0
</td>
<td style="text-align:right;">
2.7
</td>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.4
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.0
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
3.1
</td>
<td style="text-align:right;">
4.7
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.3
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:right;">
4.4
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
4.1
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.5
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:right;">
4.0
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.5
</td>
<td style="text-align:right;">
2.6
</td>
<td style="text-align:right;">
4.4
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
4.6
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.8
</td>
<td style="text-align:right;">
2.6
</td>
<td style="text-align:right;">
4.0
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:right;">
3.3
</td>
<td style="text-align:right;">
1.0
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
2.7
</td>
<td style="text-align:right;">
4.2
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
4.2
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
4.2
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.2
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
4.3
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
1.1
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
4.1
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:left;">
versicolor
</td>
</tr>
<tr>
<td style="text-align:right;">
6.3
</td>
<td style="text-align:right;">
3.3
</td>
<td style="text-align:right;">
6.0
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
5.8
</td>
<td style="text-align:right;">
2.7
</td>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
7.1
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
5.9
</td>
<td style="text-align:right;">
2.1
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.3
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.5
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
5.8
</td>
<td style="text-align:right;">
2.2
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
7.6
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
6.6
</td>
<td style="text-align:right;">
2.1
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
4.9
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:right;">
4.5
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
7.3
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
6.3
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:right;">
5.8
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
7.2
</td>
<td style="text-align:right;">
3.6
</td>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.5
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
2.0
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.4
</td>
<td style="text-align:right;">
2.7
</td>
<td style="text-align:right;">
5.3
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.8
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
5.5
</td>
<td style="text-align:right;">
2.1
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
2.0
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
5.8
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
2.4
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.4
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
5.3
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.5
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
5.5
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
7.7
</td>
<td style="text-align:right;">
3.8
</td>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
2.2
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
7.7
</td>
<td style="text-align:right;">
2.6
</td>
<td style="text-align:right;">
6.9
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.0
</td>
<td style="text-align:right;">
2.2
</td>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.9
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
4.9
</td>
<td style="text-align:right;">
2.0
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
7.7
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
2.0
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.3
</td>
<td style="text-align:right;">
2.7
</td>
<td style="text-align:right;">
4.9
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
3.3
</td>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
2.1
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
7.2
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
6.0
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.2
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
4.8
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
4.9
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.4
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
2.1
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
7.2
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
5.8
</td>
<td style="text-align:right;">
1.6
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
7.4
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
7.9
</td>
<td style="text-align:right;">
3.8
</td>
<td style="text-align:right;">
6.4
</td>
<td style="text-align:right;">
2.0
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.4
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
2.2
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.3
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
2.6
</td>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
7.7
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.3
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
2.4
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.4
</td>
<td style="text-align:right;">
3.1
</td>
<td style="text-align:right;">
5.5
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.0
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
4.8
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.9
</td>
<td style="text-align:right;">
3.1
</td>
<td style="text-align:right;">
5.4
</td>
<td style="text-align:right;">
2.1
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
3.1
</td>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
2.4
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.9
</td>
<td style="text-align:right;">
3.1
</td>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
5.8
</td>
<td style="text-align:right;">
2.7
</td>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.8
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
5.9
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
3.3
</td>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
5.2
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.3
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.5
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
5.2
</td>
<td style="text-align:right;">
2.0
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
6.2
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
5.4
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
<tr>
<td style="text-align:right;">
5.9
</td>
<td style="text-align:right;">
3.0
</td>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:left;">
virginica
</td>
</tr>
</tbody>
</table>

You can also use “inline” code, e.g. I might want to say that the iris
table has 150 rows and 5 columns, and that the mean of the Sepal.Length
column is 5.8433333

# Other useful tips

Any time you’re copy-pasting multiple lines of code, you probably should
be using a function. Well worth your time to learn how to do that!

Naming your variables - it’s very tempting to use anonymous-sounding
variable names, like `temp` or `results` or `output`. But that is likely
to cause you trouble in the long run - when you come back and look at
your code a year from now, you won’t understand it. Also, you are likely
to accidentally have two different things in the same script with the
same name - that’ll mess up your results.

# Reproducibility

If you’re really being strict about code reproducibility, it’s nice to
see which version of R and which versions of the packages were used to
generate a report. Mostly code runs in a stable way, but once in a
while, you update R or certain packages, and it either breaks the code,
or you get output that’s mysteriously different. In that case it can
help to include R and package versions in the report, like this
(typically at the end of the report). I don’t often do this but maybe I
should.

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: x86_64-apple-darwin20
    ## Running under: macOS Ventura 13.6.6
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
    ##  [1] kableExtra_1.4.0 lubridate_1.9.3  forcats_1.0.0    stringr_1.5.1   
    ##  [5] dplyr_1.1.4      purrr_1.0.2      readr_2.1.5      tidyr_1.3.1     
    ##  [9] tibble_3.2.1     ggplot2_3.5.1    tidyverse_2.0.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.5      highr_0.10        compiler_4.4.0    tidyselect_1.2.1 
    ##  [5] xml2_1.3.6        systemfonts_1.0.6 scales_1.3.0      yaml_2.3.8       
    ##  [9] fastmap_1.1.1     R6_2.5.1          generics_0.1.3    knitr_1.46       
    ## [13] munsell_0.5.1     svglite_2.1.3     pillar_1.9.0      tzdb_0.4.0       
    ## [17] rlang_1.1.3       utf8_1.2.4        stringi_1.8.4     xfun_0.43        
    ## [21] viridisLite_0.4.2 timechange_0.3.0  cli_3.6.2         withr_3.0.0      
    ## [25] magrittr_2.0.3    digest_0.6.35     grid_4.4.0        rstudioapi_0.16.0
    ## [29] hms_1.1.3         lifecycle_1.0.4   vctrs_0.6.5       evaluate_0.23    
    ## [33] glue_1.7.0        fansi_1.0.6       colorspace_2.1-0  rmarkdown_2.26   
    ## [37] tools_4.4.0       pkgconfig_2.0.3   htmltools_0.5.8.1
