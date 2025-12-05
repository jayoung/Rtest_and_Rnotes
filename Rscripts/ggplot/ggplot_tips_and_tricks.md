ggplot_tips_and_tricks
================
Janet Young

2025-12-05

# Goal

A place to collect various ggplot tips and tricks

xxxx maybe I should move the stuff in the `notes/ggplot_notes/` folder
here

``` r
## the above is a good chunk header for chunks that load libraries
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(patchwork)
library(here)
library(janitor)
library(kableExtra)

library(ggExtra)  ## for ggMarginal
library(shadowtext)
library(ggtext)
library(ggbreak) 

library(palmerpenguins)

### for my_ggMarginal:
source(here("useful_functions/other_functions.R"))
```

# Public ggplot2 resources

[ggplot2 gallery](https://www.r-graph-gallery.com/ggplot2-package.html)

[ggplot2 extensions
gallery](https://exts.ggplot2.tidyverse.org/gallery/) - shows various
packages that extent ggplot. Examples (but there are many many more): -
patchwork and cowplot to combine plots - gganimate - ggstatsplot (and
ggsignif) to add results of statistical tests, - gggenomes and gggenes

A ggplot
[tutorial](https://rgup.gitlab.io/research_cycle/02_ggplot.html)

[Demos](https://www.datanovia.com/en/lessons/combine-multiple-ggplots-into-a-figure/#google_vignette)
of how to combine \>1 plot into a figure. Packages:

- gridExtra::grid.arrange()
- cowplot::plot_grid()
- patchwork::plot_layout()
- ggpubr::ggarrange()

<https://www.thoughtworks.com/insights/blog/coding-habits-data-scientists>

# ggplot2 version notes

ggplot2 version 4.0.0 is a big change, and may break some things in
other packages. See notes from the Bioconductor team
[here](https://blog.bioconductor.org/posts/2025-07-07-ggplot2-update/).

I don’t want to update to it yet. That may mean I need to use older
versions of some other packages.

Here’s how I can install a particular archive version of a package once
I identify the URL on CRAN:

``` r
# ### revert back to ggplot 3.5.2 - this is from 2025-04-09 
# ## version 4.0.0 is
# packageurl <- "http://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.5.2.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")

### revert back to GGally 2.2.1 (2024-02-14). Newer versions need ggplot2 version 4.4.0
# # GGally_2.2.0.tar.gz
# packageurl <- "http://cran.r-project.org/src/contrib/Archive/GGally/GGally_2.2.1.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
```

# ggplot tips and tricks

## Example datasets for this script

Get a couple of example datasets - show the first few rows of each

### iris_tbl

(clean up the built-in iris dataset a bit)

``` r
iris_tbl <- iris %>% 
    as_tibble() %>% 
    clean_names()
```

`iris_tbl` has 150 rows and 5 columns

``` r
iris_tbl %>% 
    slice_head(n=3) %>% 
    kable(caption="iris_tbl") %>% 
    kable_styling(full_width = FALSE)
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

iris_tbl
</caption>

<thead>

<tr>

<th style="text-align:right;">

sepal_length
</th>

<th style="text-align:right;">

sepal_width
</th>

<th style="text-align:right;">

petal_length
</th>

<th style="text-align:right;">

petal_width
</th>

<th style="text-align:left;">

species
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

</tbody>

</table>

### penguins

`penguins` has 344 rows and 8 columns

``` r
penguins %>% 
    slice_head(n=3) %>% 
    kable(caption="penguins") %>% 
    kable_styling(full_width = FALSE)
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

penguins
</caption>

<thead>

<tr>

<th style="text-align:left;">

species
</th>

<th style="text-align:left;">

island
</th>

<th style="text-align:right;">

bill_length_mm
</th>

<th style="text-align:right;">

bill_depth_mm
</th>

<th style="text-align:right;">

flipper_length_mm
</th>

<th style="text-align:right;">

body_mass_g
</th>

<th style="text-align:left;">

sex
</th>

<th style="text-align:right;">

year
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Adelie
</td>

<td style="text-align:left;">

Torgersen
</td>

<td style="text-align:right;">

39.1
</td>

<td style="text-align:right;">

18.7
</td>

<td style="text-align:right;">

181
</td>

<td style="text-align:right;">

3750
</td>

<td style="text-align:left;">

male
</td>

<td style="text-align:right;">

2007
</td>

</tr>

<tr>

<td style="text-align:left;">

Adelie
</td>

<td style="text-align:left;">

Torgersen
</td>

<td style="text-align:right;">

39.5
</td>

<td style="text-align:right;">

17.4
</td>

<td style="text-align:right;">

186
</td>

<td style="text-align:right;">

3800
</td>

<td style="text-align:left;">

female
</td>

<td style="text-align:right;">

2007
</td>

</tr>

<tr>

<td style="text-align:left;">

Adelie
</td>

<td style="text-align:left;">

Torgersen
</td>

<td style="text-align:right;">

40.3
</td>

<td style="text-align:right;">

18.0
</td>

<td style="text-align:right;">

195
</td>

<td style="text-align:right;">

3250
</td>

<td style="text-align:left;">

female
</td>

<td style="text-align:right;">

2007
</td>

</tr>

</tbody>

</table>

## Using color names in a column literally

Use `scale_color_identity()` if colors are specified by name in a
column. Note that by default `guide="none"` for `scale_color_identity()`
(could add `guide=guide_legend()` but that would just show the color
names, which makes no sense)

``` r
p1 <- iris_tbl %>% 
    ggplot(aes(x=sepal_length, y=petal_length, color=species)) +
    geom_point() + 
    theme_classic() +
    labs(title="Color points by species")

p2 <- iris_tbl %>% 
    mutate(species_color=case_when(
        species=="versicolor" ~ "orange",
        TRUE ~ "darkgray"
    )) %>% 
    ggplot(aes(x=sepal_length, y=petal_length, color=species_color)) +
    geom_point() + 
    theme_classic() +
    scale_color_identity() +
    labs(title="Color to highlight versicolor")

p1 + p2
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Adding individual data points to grouped/filled boxplots and getting the aligned

The trick is that we need to ‘dodge’ the points, using
`geom_point(position=position_jitterdodge())` (or we could have used
`position_dodge` if we don’t want the jitter)

``` r
iris_plus_groups <- iris_tbl %>% 
    as_tibble() %>% 
    mutate(group = sample(1:2, size=nrow(iris_tbl), replace=TRUE)) %>% 
    mutate(group=paste0("group_", group)) 

iris_plus_groups %>% 
    ggplot(aes(x=species, y=sepal_length, color=group)) +
    geom_boxplot() +
    geom_point(position=position_jitterdodge(jitter.width = 0.12)) +
    theme_classic()
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

If there’s an empty group, widths and spacing get weird (see left plot
below), so we do
`geom_boxplot(position = position_dodge(preserve = "single"))` (see
right plot below).

``` r
p1 <- iris_plus_groups %>% 
    filter( ! (species=="versicolor" & group=="group_2")  ) %>% 
    ggplot(aes(x=species, y=sepal_length, color=group)) +
    geom_boxplot() +
    theme_classic()
p2 <- iris_plus_groups %>% 
    filter( ! (species=="versicolor" & group=="group_2")  ) %>% 
    ggplot(aes(x=species, y=sepal_length, color=group)) +
    geom_boxplot(position = position_dodge(preserve = "single")) +
    theme_classic()
(p1 + p2) +
    plot_layout(guides="collect")
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

The `position_dodge2` function gives ALMOST the same outputas
`position_dodge`, but the alignment of the ‘versicolor’ box is
different:

``` r
iris_plus_groups %>% 
    filter( ! (species=="versicolor" & group=="group_2")  ) %>% 
    ggplot(aes(x=species, y=sepal_length, color=group)) +
    geom_boxplot(position = position_dodge2(preserve = "single")) +
    theme_classic()
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

To ALSO add geom_point and keep them lined up is tricky! This might be
fixed in newer versions of ggplot2 - there are a few related bug
reports. In the meantime this is a workaround (given
[here](https://github.com/tidyverse/ggplot2/issues/2712)) - we have to
use `position_dodge2` for the boxplots and `position_jitterdodge` for
the points:

``` r
iris_plus_groups %>% 
    filter( ! (species=="versicolor" & group=="group_2")  ) %>% 
    ggplot(aes(x=species, y=sepal_length, color=group)) + 
    geom_boxplot(position = position_dodge2(0.75, preserve = 'single')) +
    geom_point(position = position_jitterdodge(dodge.width=0.75, jitter.width=0.15))+
    theme_classic()
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Violin plots

### Things I usually do on violin plots

Add median dots:

    + stat_summary(fun = "median", geom = "point",
                   position=position_dodge(width=0.9),
                   show.legend = FALSE)  

Rotate x axis labels:

`+ theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))`

``` r
penguins %>%
    filter(!is.na(flipper_length_mm)) %>% 
    ggplot(aes(x=species, y=flipper_length_mm, fill=species)) +
    geom_violin(scale = "width") + 
    stat_summary(fun = "median", geom = "point",
                 ## position_dodge(width=0.9) will help median dots line up especially if each x axis value contains subgroups
                 position=position_dodge(width=0.9),
                 show.legend = FALSE)  +
    theme_classic() +
    labs(x="") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    guides(fill="none")
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Note: I sometimes have trouble getting the statistical summary (median
dot) to line up correctly if I’ve got subgroupings for the violins. The
trick to put `position=position_dodge(width=0.9)` in the stat_summary
call.

## Annotating plots - `annotate()` versus `geom_text()`

`annotate()` is better than `geom_text()` for some uses.

Demo based on
[rfortherestofus](https://rfortherestofus.com/2023/10/annotate-vs-geoms)

- use geom_text if the data itself should drive text labels
- use annotate if you’re manually adding labels

Example where annotate is better

``` r
scatterplot <- palmerpenguins::penguins %>% 
    select(bill_length_mm, flipper_length_mm, species) %>% 
    drop_na()%>% 
    ggplot(aes(bill_length_mm, flipper_length_mm, col = species)) +
    geom_point(size = 2.5) +
    labs(
        x = 'Bill length (in mm)',
        y = 'Flipper length (in mm)',
        col = 'species',
        title = 'Measurements of different Penguin species'
    ) +
    theme_minimal(base_size = 16) +
    theme(legend.position = 'top')
```

Here we use `geom_text()` and it looks bad, because it is still drawing
stuff (color) from the data passed in, and that’s not what we want

``` r
scatterplot +
    geom_text(
        x = 35,
        y = 217.5,
        label = 'Important penguins',
        fontface = 'bold', # makes text bold
        size = 4.5 # font size
    )
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Here we use `annotate()` instead and it ignores the data

``` r
scatterplot +
    annotate(
        'text',
        x = 35,
        y = 217.5,
        label = 'Important penguins',
        fontface = 'bold', 
        size = 4.5
    )
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

## Explore `shadowtext` package

Demo from
[rfortherestofus](https://rfortherestofus.com/2024/05/shadowtext-ggplot)

Mostly I don’t like the way shadowtext looks, but it could be really
useful if we’re putting words over the top of some other dense data,
e.g. words on top of a map benefit from a white shadow

``` r
### revert back to shadowtext 0.1.4 (2024-07-18). Newer versions need ggplot2 version 4.4.0
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/shadowtext/shadowtext_0.1.4.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
```

``` r
### to install missing fonts:
# library(showtext) # For fonts
# font_add_google("Source Sans Pro") ## downloads and installs a font, except I don't think it necessarily makes the font available always?  maybe I need to restart the computer before I can use it?
# grep("Source", system_fonts()$name, value=TRUE)
# grep("Times", system_fonts()$name, value=TRUE)
# grep("Arial", system_fonts()$name, value=TRUE)
```

``` r
## define labelling info
species_labels_tib <- tibble(
    species = c('Adelie', 'Gentoo', 'Chinstrap'),
    x = c(35, 43, 53),
    y = c(210, 229, 178)
)
```

We can use geom_shadowtext() pretty much the same as we’d use
geom_label, using bg.color to tell it the color of the shadow.

``` r
## the plot
penguins %>% 
    drop_na() %>%
    ggplot(
        aes(bill_length_mm, flipper_length_mm, fill = species)
    ) +
    geom_point(shape = 21, size = 2) +
    geom_shadowtext(
        data = species_labels_tib,
        aes(x, y, col = species, label = species),
        size = 6,
        fontface = 'bold',
        family = 'ArialMT',
        bg.color = 'grey10',
    )  +
    theme_minimal(base_size = 16, base_family = 'ArialMT') +
    theme(legend.position = 'none')
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

## Wrapping text over \>1 line in ggplot

`str_wrap()` - you have to figure out width manually, which can be
tedious

``` r
penguins %>% 
    count(island) %>%
    ggplot(aes(x=island, y=n)) +
    geom_col() +
    labs(title="a short title",
         subtitle=str_wrap("a really long title. kasjdhf ;isjdghf khg kajsxdhf khg alsidgf kjhg ljhags dfj hgkjahsdgfkjhg a  ljhsdgf ljhglsdjhfg", width=50))
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Instead use `ggtext` package - the element_textbox_simple will
automatically wrap text to fit whatever space is available.

``` r
penguins %>% 
    count(island) %>%
    ggplot(aes(x=island, y=n)) +
    geom_col() +
    labs(title="a short title",
         subtitle="a really long title. kasjdhf ;isjdghf khg kajsdhf khg alsidgf kjhg ljhags dfj hgkjahsdgfkjhg a  ljhsdgf ljhglsdjhfg") +
    theme(plot.subtitle = element_textbox_simple())
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

Maybe we need to wrap facet labels - we can use the label_wrap_gen
function

``` r
theme_set(theme_bw(16))
df <- data.frame(measurement = rnorm(20,mean=30), 
                 group = c(rep('A really long group variable name that needs to be wrapped',10),
                           rep('group1',10)), 
                 sex = c(rep(c('M','F'),10)))

df %>%
    ggplot(aes(sex, measurement, color = sex)) +
    geom_boxplot() +
    facet_wrap(~group, labeller = label_wrap_gen(width=24))+
    theme(legend.position="none")
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

## Discontinuous axes using `ggbreak` package

`ggbreak` package
[vignette](https://cran.r-project.org/web/packages/ggbreak/vignettes/ggbreak.html)

there is a blank plot below, as well as the intended plots, but it
doesn’t appear when you knit to html or github_document

``` r
set.seed(2019-01-19)
d <- data.frame(x = 1:20,
                y = c(rnorm(5) + 4, rnorm(5) + 20, rnorm(5) + 5, rnorm(5) + 22)
)

p1 <- ggplot(d, aes(y, x)) + 
    geom_col(orientation="y") +
    theme_minimal() +
    labs(title="ordinary x-axis")
p2 <- p1 + 
    scale_x_break(c(7, 17)) + # from ggbreak
    labs(title="broken x-axis")


p1 + p2
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

## Add marginal density plots

We use [`ggExtra::ggMarginal`](https://github.com/daattali/ggExtra) to
add marginal density plots.

We do it a slightly weird way because it’s a ggExtraPlot object not a
regular ggplot object - see [this
note](https://github.com/daattali/ggExtra?tab=readme-ov-file#using-ggmarginal-in-r-notebooks-or-rmarkdown)

``` r
## first make a basic plot
p1 <- penguins %>% 
    drop_na() %>%
    ggplot(
        aes(bill_length_mm, flipper_length_mm, color = species)
    ) +
    geom_point(size = 1) +
    theme_classic() + 
    ## legend needs to be at bottom (or left) so it doesn't push the density plot away from the main plot
    theme(legend.position="bottom") 

## then use ggMarginal
p1a <- ggMarginal(p1, type="density", groupColour = TRUE)
```

``` r
p1a
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

If we want to use patchwork to combine \>1 ggMarginal plot, we need to
use the `wrap_elements()` function:

``` r
patchwork::wrap_elements(p1a) + patchwork::wrap_elements(p1a)
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

[A known issue with
ggMarginal](https://github.com/daattali/ggExtra/issues/128) - it doesn’t
respect coord_cartesian settings on the main plot:

``` r
p1b <- p1 +
    coord_cartesian(xlim=c(0,60)) 
p1b <- ggMarginal(p1b, type="density", groupColour = TRUE)
```

``` r
p1b
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

A workaround is to use `xlim` rather than coord_cartesian. The downside
of that is that it actually REMOVES data outside the specified range, so
the density plots are not accurate.

``` r
p1b <- p1 +
    xlim(c(0,60))
p1b <- ggMarginal(p1b, type="density", groupColour = TRUE)
```

``` r
p1b
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

We can also do it in basic ggplot (including changing the axis limits).
I made a function to do that - it’s called `my_ggMarginal()` and is in
`useful_functions/other_functions.R`

Show that function:

``` r
print(my_ggMarginal)
```

    ## function (df, x_var = NULL, y_var = NULL, color_var = NULL, my_xlim = NULL, 
    ##     my_ylim = NULL, my_color_scheme = NULL, my_title = NULL, 
    ##     my_subtitle = NULL, combine_plots = TRUE, ...) 
    ## {
    ##     if (is.null(x_var) | is.null(y_var)) {
    ##         stop("\n\nYou must supply x_var and y_var\n\n")
    ##     }
    ##     if (!x_var %in% colnames(df)) {
    ##         stop("\n\nYour x variable name ", x_var, " is not in the data frame you supplied\n\n")
    ##     }
    ##     if (!y_var %in% colnames(df)) {
    ##         stop("\n\nYour y variable name ", y_var, " is not in the data frame you supplied\n\n")
    ##     }
    ##     if (!color_var %in% colnames(df)) {
    ##         stop("\n\nYour colors variable name ", color_var, " is not in the data frame you supplied\n\n")
    ##     }
    ##     if (is.null(color_var)) {
    ##         p1 <- df %>% ggplot(aes(x = .data[[x_var]], y = .data[[y_var]]))
    ##         dens_x <- df %>% ggplot(aes(x = .data[[x_var]]))
    ##         dens_y <- df %>% ggplot(aes(y = .data[[y_var]]))
    ##     }
    ##     else {
    ##         p1 <- df %>% ggplot(aes(x = .data[[x_var]], y = .data[[y_var]], 
    ##             color = .data[[color_var]]))
    ##         dens_x <- df %>% ggplot(aes(x = .data[[x_var]], color = .data[[color_var]], 
    ##             fill = .data[[color_var]]))
    ##         dens_y <- df %>% ggplot(aes(y = .data[[y_var]], color = .data[[color_var]], 
    ##             fill = .data[[color_var]]))
    ##     }
    ##     p1 <- p1 + geom_point(...) + theme_classic()
    ##     dens_x <- dens_x + geom_density(show.legend = FALSE, alpha = 0.3) + 
    ##         theme_void()
    ##     dens_y <- dens_y + geom_density(show.legend = FALSE, alpha = 0.3) + 
    ##         theme_void()
    ##     if (!is.null(my_xlim) & !is.null(my_ylim)) {
    ##         p1 <- p1 + coord_cartesian(xlim = my_xlim, ylim = my_ylim)
    ##     }
    ##     if (!is.null(my_xlim)) {
    ##         if (is.null(my_ylim)) {
    ##             p1 <- p1 + coord_cartesian(xlim = my_xlim)
    ##         }
    ##         dens_x <- dens_x + coord_cartesian(xlim = my_xlim)
    ##     }
    ##     if (!is.null(my_ylim)) {
    ##         if (is.null(my_xlim)) {
    ##             p1 <- p1 + coord_cartesian(ylim = my_ylim)
    ##         }
    ##         dens_y <- dens_y + coord_cartesian(ylim = my_ylim)
    ##     }
    ##     if (!is.null(my_color_scheme)) {
    ##         p1 <- p1 + scale_color_manual(values = my_color_scheme)
    ##         dens_x <- dens_x + scale_color_manual(values = my_color_scheme)
    ##         dens_y <- dens_y + scale_color_manual(values = my_color_scheme)
    ##     }
    ##     if (!is.null(my_title)) {
    ##         dens_x <- dens_x + labs(title = my_title)
    ##     }
    ##     if (!is.null(my_subtitle)) {
    ##         dens_x <- dens_x + labs(subtitle = my_subtitle)
    ##     }
    ##     if (!combine_plots) {
    ##         return(list(main_plot = p1, dens_x = dens_x, dens_y = dens_y))
    ##     }
    ##     all_plots <- (dens_x + plot_spacer() + p1 + dens_y) + plot_layout(ncol = 2, 
    ##         nrow = 2, widths = c(6, 1), heights = c(1, 6), guides = "collect")
    ##     return(all_plots)
    ## }

``` r
my_penguin_colors <- c(Adelie="orange",
                       Gentoo="purple",
                       Chinstrap="forestgreen")
penguins %>%  
    drop_na() %>% 
    my_ggMarginal(x_var="bill_length_mm",
                  y_var="flipper_length_mm", color_var="species",
                  my_xlim=c(0,60), my_ylim=c(160,240),
                  my_color_scheme=my_penguin_colors,
                  my_title="penguins",
                  my_subtitle = "by species")
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

If we want to show \>1 of those plots together (from `my_ggMarginal`),
we need a special function. First we add the `combine_plots = FALSE`
option when we run `my_ggMarginal()`.

Then we combine plots into a list, and feed that to
`combine_myGGmarginal_plots`:

``` r
p1 <- penguins %>%  
    drop_na() %>% 
    my_ggMarginal(x_var="bill_length_mm",
                  y_var="flipper_length_mm", 
                  color_var="species",
                  my_xlim=c(0,60), my_ylim=c(160,240),
                  my_color_scheme=my_penguin_colors,
                  my_title="penguins",
                  my_subtitle = "by species", 
                  combine_plots = FALSE)

combine_myGGmarginal_plots(list(p1,p1))
```

![](ggplot_tips_and_tricks_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

`geom_rug()` is an alternative way to show the marginal distributions:

``` r
p1 +
    geom_rug(alpha=.2,
             length = unit(0.1, "inches")) 
```

    ## NULL

``` r
## default length is unit(0.03, "npc"), i.e. 0.03* the plot dimensions (so the x and y axis rugs might be different lengths, unless we control that)
```

# Finished

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Tahoe 26.1
    ## 
    ## Matrix products: default
    ## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
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
    ##  [1] palmerpenguins_0.1.1 ggbreak_0.1.6        ggtext_0.1.2        
    ##  [4] shadowtext_0.1.4     ggExtra_0.11.0       kableExtra_1.4.0    
    ##  [7] janitor_2.2.1        here_1.0.2           patchwork_1.3.2     
    ## [10] lubridate_1.9.4      forcats_1.0.0        stringr_1.5.2       
    ## [13] dplyr_1.1.4          purrr_1.1.0          readr_2.1.5         
    ## [16] tidyr_1.3.1          tibble_3.3.0         ggplot2_3.5.2       
    ## [19] tidyverse_2.0.0     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6       xfun_0.53          tzdb_0.5.0         vctrs_0.6.5       
    ##  [5] tools_4.5.2        generics_0.1.4     yulab.utils_0.2.1  pkgconfig_2.0.3   
    ##  [9] ggplotify_0.1.3    RColorBrewer_1.1-3 lifecycle_1.0.4    compiler_4.5.2    
    ## [13] farver_2.1.2       textshaping_1.0.3  snakecase_0.11.1   litedown_0.7      
    ## [17] ggfun_0.2.0        httpuv_1.6.16      htmltools_0.5.8.1  yaml_2.3.10       
    ## [21] pillar_1.11.1      later_1.4.4        mime_0.13          commonmark_2.0.0  
    ## [25] aplot_0.2.9        tidyselect_1.2.1   digest_0.6.37      stringi_1.8.7     
    ## [29] labeling_0.4.3     rprojroot_2.1.1    fastmap_1.2.0      grid_4.5.2        
    ## [33] cli_3.6.5          magrittr_2.0.4     withr_3.0.2        scales_1.4.0      
    ## [37] promises_1.3.3     rappdirs_0.3.3     timechange_0.3.0   rmarkdown_2.29    
    ## [41] hms_1.1.3          shiny_1.11.1       evaluate_1.0.5     knitr_1.50        
    ## [45] miniUI_0.1.2       viridisLite_0.4.2  markdown_2.0       gridGraphics_0.5-1
    ## [49] rlang_1.1.6        gridtext_0.1.5     Rcpp_1.1.0         xtable_1.8-4      
    ## [53] glue_1.8.0         xml2_1.4.0         svglite_2.2.1      rstudioapi_0.17.1 
    ## [57] R6_2.6.1           systemfonts_1.3.1  fs_1.6.6
