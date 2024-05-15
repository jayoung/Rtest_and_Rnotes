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
iris
```

    ##     Sepal.Length Sepal.Width Petal.Length Petal.Width    Species
    ## 1            5.1         3.5          1.4         0.2     setosa
    ## 2            4.9         3.0          1.4         0.2     setosa
    ## 3            4.7         3.2          1.3         0.2     setosa
    ## 4            4.6         3.1          1.5         0.2     setosa
    ## 5            5.0         3.6          1.4         0.2     setosa
    ## 6            5.4         3.9          1.7         0.4     setosa
    ## 7            4.6         3.4          1.4         0.3     setosa
    ## 8            5.0         3.4          1.5         0.2     setosa
    ## 9            4.4         2.9          1.4         0.2     setosa
    ## 10           4.9         3.1          1.5         0.1     setosa
    ## 11           5.4         3.7          1.5         0.2     setosa
    ## 12           4.8         3.4          1.6         0.2     setosa
    ## 13           4.8         3.0          1.4         0.1     setosa
    ## 14           4.3         3.0          1.1         0.1     setosa
    ## 15           5.8         4.0          1.2         0.2     setosa
    ## 16           5.7         4.4          1.5         0.4     setosa
    ## 17           5.4         3.9          1.3         0.4     setosa
    ## 18           5.1         3.5          1.4         0.3     setosa
    ## 19           5.7         3.8          1.7         0.3     setosa
    ## 20           5.1         3.8          1.5         0.3     setosa
    ## 21           5.4         3.4          1.7         0.2     setosa
    ## 22           5.1         3.7          1.5         0.4     setosa
    ## 23           4.6         3.6          1.0         0.2     setosa
    ## 24           5.1         3.3          1.7         0.5     setosa
    ## 25           4.8         3.4          1.9         0.2     setosa
    ## 26           5.0         3.0          1.6         0.2     setosa
    ## 27           5.0         3.4          1.6         0.4     setosa
    ## 28           5.2         3.5          1.5         0.2     setosa
    ## 29           5.2         3.4          1.4         0.2     setosa
    ## 30           4.7         3.2          1.6         0.2     setosa
    ## 31           4.8         3.1          1.6         0.2     setosa
    ## 32           5.4         3.4          1.5         0.4     setosa
    ## 33           5.2         4.1          1.5         0.1     setosa
    ## 34           5.5         4.2          1.4         0.2     setosa
    ## 35           4.9         3.1          1.5         0.2     setosa
    ## 36           5.0         3.2          1.2         0.2     setosa
    ## 37           5.5         3.5          1.3         0.2     setosa
    ## 38           4.9         3.6          1.4         0.1     setosa
    ## 39           4.4         3.0          1.3         0.2     setosa
    ## 40           5.1         3.4          1.5         0.2     setosa
    ## 41           5.0         3.5          1.3         0.3     setosa
    ## 42           4.5         2.3          1.3         0.3     setosa
    ## 43           4.4         3.2          1.3         0.2     setosa
    ## 44           5.0         3.5          1.6         0.6     setosa
    ## 45           5.1         3.8          1.9         0.4     setosa
    ## 46           4.8         3.0          1.4         0.3     setosa
    ## 47           5.1         3.8          1.6         0.2     setosa
    ## 48           4.6         3.2          1.4         0.2     setosa
    ## 49           5.3         3.7          1.5         0.2     setosa
    ## 50           5.0         3.3          1.4         0.2     setosa
    ## 51           7.0         3.2          4.7         1.4 versicolor
    ## 52           6.4         3.2          4.5         1.5 versicolor
    ## 53           6.9         3.1          4.9         1.5 versicolor
    ## 54           5.5         2.3          4.0         1.3 versicolor
    ## 55           6.5         2.8          4.6         1.5 versicolor
    ## 56           5.7         2.8          4.5         1.3 versicolor
    ## 57           6.3         3.3          4.7         1.6 versicolor
    ## 58           4.9         2.4          3.3         1.0 versicolor
    ## 59           6.6         2.9          4.6         1.3 versicolor
    ## 60           5.2         2.7          3.9         1.4 versicolor
    ## 61           5.0         2.0          3.5         1.0 versicolor
    ## 62           5.9         3.0          4.2         1.5 versicolor
    ## 63           6.0         2.2          4.0         1.0 versicolor
    ## 64           6.1         2.9          4.7         1.4 versicolor
    ## 65           5.6         2.9          3.6         1.3 versicolor
    ## 66           6.7         3.1          4.4         1.4 versicolor
    ## 67           5.6         3.0          4.5         1.5 versicolor
    ## 68           5.8         2.7          4.1         1.0 versicolor
    ## 69           6.2         2.2          4.5         1.5 versicolor
    ## 70           5.6         2.5          3.9         1.1 versicolor
    ## 71           5.9         3.2          4.8         1.8 versicolor
    ## 72           6.1         2.8          4.0         1.3 versicolor
    ## 73           6.3         2.5          4.9         1.5 versicolor
    ## 74           6.1         2.8          4.7         1.2 versicolor
    ## 75           6.4         2.9          4.3         1.3 versicolor
    ## 76           6.6         3.0          4.4         1.4 versicolor
    ## 77           6.8         2.8          4.8         1.4 versicolor
    ## 78           6.7         3.0          5.0         1.7 versicolor
    ## 79           6.0         2.9          4.5         1.5 versicolor
    ## 80           5.7         2.6          3.5         1.0 versicolor
    ## 81           5.5         2.4          3.8         1.1 versicolor
    ## 82           5.5         2.4          3.7         1.0 versicolor
    ## 83           5.8         2.7          3.9         1.2 versicolor
    ## 84           6.0         2.7          5.1         1.6 versicolor
    ## 85           5.4         3.0          4.5         1.5 versicolor
    ## 86           6.0         3.4          4.5         1.6 versicolor
    ## 87           6.7         3.1          4.7         1.5 versicolor
    ## 88           6.3         2.3          4.4         1.3 versicolor
    ## 89           5.6         3.0          4.1         1.3 versicolor
    ## 90           5.5         2.5          4.0         1.3 versicolor
    ## 91           5.5         2.6          4.4         1.2 versicolor
    ## 92           6.1         3.0          4.6         1.4 versicolor
    ## 93           5.8         2.6          4.0         1.2 versicolor
    ## 94           5.0         2.3          3.3         1.0 versicolor
    ## 95           5.6         2.7          4.2         1.3 versicolor
    ## 96           5.7         3.0          4.2         1.2 versicolor
    ## 97           5.7         2.9          4.2         1.3 versicolor
    ## 98           6.2         2.9          4.3         1.3 versicolor
    ## 99           5.1         2.5          3.0         1.1 versicolor
    ## 100          5.7         2.8          4.1         1.3 versicolor
    ## 101          6.3         3.3          6.0         2.5  virginica
    ## 102          5.8         2.7          5.1         1.9  virginica
    ## 103          7.1         3.0          5.9         2.1  virginica
    ## 104          6.3         2.9          5.6         1.8  virginica
    ## 105          6.5         3.0          5.8         2.2  virginica
    ## 106          7.6         3.0          6.6         2.1  virginica
    ## 107          4.9         2.5          4.5         1.7  virginica
    ## 108          7.3         2.9          6.3         1.8  virginica
    ## 109          6.7         2.5          5.8         1.8  virginica
    ## 110          7.2         3.6          6.1         2.5  virginica
    ## 111          6.5         3.2          5.1         2.0  virginica
    ## 112          6.4         2.7          5.3         1.9  virginica
    ## 113          6.8         3.0          5.5         2.1  virginica
    ## 114          5.7         2.5          5.0         2.0  virginica
    ## 115          5.8         2.8          5.1         2.4  virginica
    ## 116          6.4         3.2          5.3         2.3  virginica
    ## 117          6.5         3.0          5.5         1.8  virginica
    ## 118          7.7         3.8          6.7         2.2  virginica
    ## 119          7.7         2.6          6.9         2.3  virginica
    ## 120          6.0         2.2          5.0         1.5  virginica
    ## 121          6.9         3.2          5.7         2.3  virginica
    ## 122          5.6         2.8          4.9         2.0  virginica
    ## 123          7.7         2.8          6.7         2.0  virginica
    ## 124          6.3         2.7          4.9         1.8  virginica
    ## 125          6.7         3.3          5.7         2.1  virginica
    ## 126          7.2         3.2          6.0         1.8  virginica
    ## 127          6.2         2.8          4.8         1.8  virginica
    ## 128          6.1         3.0          4.9         1.8  virginica
    ## 129          6.4         2.8          5.6         2.1  virginica
    ## 130          7.2         3.0          5.8         1.6  virginica
    ## 131          7.4         2.8          6.1         1.9  virginica
    ## 132          7.9         3.8          6.4         2.0  virginica
    ## 133          6.4         2.8          5.6         2.2  virginica
    ## 134          6.3         2.8          5.1         1.5  virginica
    ## 135          6.1         2.6          5.6         1.4  virginica
    ## 136          7.7         3.0          6.1         2.3  virginica
    ## 137          6.3         3.4          5.6         2.4  virginica
    ## 138          6.4         3.1          5.5         1.8  virginica
    ## 139          6.0         3.0          4.8         1.8  virginica
    ## 140          6.9         3.1          5.4         2.1  virginica
    ## 141          6.7         3.1          5.6         2.4  virginica
    ## 142          6.9         3.1          5.1         2.3  virginica
    ## 143          5.8         2.7          5.1         1.9  virginica
    ## 144          6.8         3.2          5.9         2.3  virginica
    ## 145          6.7         3.3          5.7         2.5  virginica
    ## 146          6.7         3.0          5.2         2.3  virginica
    ## 147          6.3         2.5          5.0         1.9  virginica
    ## 148          6.5         3.0          5.2         2.0  virginica
    ## 149          6.2         3.4          5.4         2.3  virginica
    ## 150          5.9         3.0          5.1         1.8  virginica

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
