ggplot_demo
================
Janet Young

2024-12-18

Show some examples from the [ggplot2
cheatsheet](https://rstudio.github.io/cheatsheets/html/data-visualization.html)

These commands set up the base plots `a` and `b`, but they have no geom
layer, so nothing is actually plotted

The `economics` dataset is a 574-row 6-col tibble (wide format)
describing how 5 metrics change over time. Column names are date, pce,
pop, psavert, uempmed, unemploy

The `seals` dataset is a 1155-row 4-col tibble that describes how some
animals move. Four columns named lat, long, delta_long, delta_lat

``` r
# date, pce, pop, psavert, uempmed, unemploy
a <- ggplot(economics, aes(date,unemploy))

# lat  long delta_long delta_lat
b <- ggplot(seals, aes(x=long, y=lat))
# ?seals 
# "This vector field was produced from the data described in..."
### each longitude value in seals is present 21 times, each latitude 55 times
# seals %>% count(long) %>% count(n)
# seals %>% count(lat)  %>% count(n)
# it's simply a grid in x-y space
```

I wanted to explore different geoms

geom_blank()

``` r
a + geom_blank() + labs(title="economics, date/unemploy, geom_blank")
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
geom_point()

``` r
b + geom_point() + labs(title="seals, long/lat, geom_point")
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
geom_curve() draws a curved line from (x, y) to (xend, yend).

aes() arguments: x, xend, y, yend, alpha, angle, color, curvature,
linetype, size.

``` r
## this is a wierd plot
# x and y already existed in b, and we're adding xend and yend
b + geom_curve(aes(xend = long + 1, 
                   yend = lat + 1), 
               curvature = 1) + 
    labs(title="seals, long/lat, geom_curve")
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
`geom_path()` connects observations in the order they appear.

aes() arguments: x, y, alpha, color, group, linetype, size.

``` r
a + 
    geom_path(lineend = "butt", linejoin = "round", linemitre = 1) +
    labs(title="geom_path")
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

`geom_line()` connects points from left to right. `geom_path()` is
similar but connects points in the order they appear in the data. If
data are sorted, geom_line() and geom_path plots look the same,
otherwise they’re different.

``` r
a + 
    geom_path(lineend = "butt", linejoin = "round", linemitre = 1) +
    coord_cartesian(xlim=c(date("1970-01-01"),date("1975-01-01"))) +
    labs(title="geom_path")
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# Connect observations in the order they appear. aes() arguments: x, y, alpha, color, group, linetype, size.
```

Same thing with geom_line - looks the same (because data are sorted,
geom_path IS connected left-to-right, just like geom_line)

geom_line() - connect observations in the order they appear. a

geom_line aes() arguments: x, y, alpha, color, group, linetype, size.

``` r
a + 
    geom_line(lineend = "butt", linejoin = "round", linemitre = 1) +
    coord_cartesian(xlim=c(date("1970-01-01"),date("1975-01-01"))) +
    labs(title="geom_line")
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

`geom_polygon()` - connect points into polygons

aes() arguments: x, y, alpha, color, fill, group, subgroup, linetype,
size.

``` r
a + geom_polygon(aes(alpha = 50)) +
    labs(title="geom_polygon")
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

`geom_rect()` draws a rectangle by connecting four corners (xmin, xmax,
ymin, ymax).

aes() arguments: xmax, xmin, ymax, ymin, alpha, color, fill, linetype,
size.

the example given makes it hard to see what’s going on, so I sample 30
rows at random

``` r
seals %>% 
    slice_sample(n=30) %>% 
    ggplot(aes(xmin = long, ymin = lat, 
               xmax = long + 1, ymax = lat + 1)) + 
    geom_rect() +
    labs(title="geom_rect")
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

`geom_ribbon()` - for each x, plot an interval from ymin to ymax.

aes() arguments: x, ymax, ymin, alpha, color, fill, group, linetype,
size.

``` r
a + geom_ribbon(aes(ymin = unemploy - 900, 
                    ymax = unemploy + 900)) +
    labs(title="geom_ribbon")
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-10-1.png)<!-- --> Plots
using `mpg` dataset - a 234 × 11 tibble, showing data for various car
types

``` r
e <- ggplot(mpg, aes(cty,hwy))
```

`geom_point()`

``` r
e + 
    geom_point() +
    labs(title="mpg, geom_point")
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

geom_smooth()

``` r
e + 
    geom_point() + 
    geom_smooth()+
    labs(title="geom_point and geom_smooth")
```

    ## `geom_smooth()` using method = 'loess' and formula = 'y ~ x'

![](ggplot_demo_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

# shadowtext package

Demo from
[rfortherestofus](https://rfortherestofus.com/2024/05/shadowtext-ggplot)

Mostly I don’t like the way shadowtext looks, but it could be really
useful if we’re putting words over the top of some other dense data,
e.g. words on top of a map benefit from a white shadow

``` r
### to install missing fonts:
# library(showtext) # For fonts
# font_add_google("Source Sans Pro") ## downloads and installs a font, except I don't think it necessarily makes the font available always?  maybe I need to restart the computer before I can use it?
# grep("Source", system_fonts()$name, value=TRUE)
# grep("Times", system_fonts()$name, value=TRUE)
# grep("Arial", system_fonts()$name, value=TRUE)
```

``` r
library(shadowtext)
library(palmerpenguins)
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
penguins |> 
    drop_na() |>
  ggplot(
    aes(bill_length_mm, flipper_length_mm, fill = species)
  ) +
  geom_point(shape = 21, size = 4) +
  geom_shadowtext(
    data = species_labels_tib,
    aes(x, y, col = species, label = species),
    size = 12,
    fontface = 'bold',
    family = 'ArialMT',
    bg.color = 'grey10',
  )  +
  theme_minimal(base_size = 16, base_family = 'ArialMT') +
  theme(legend.position = 'none')
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

# wrapping text in ggplot

`str_wrap()` - you have to figure out width manually, which can be
tedious

``` r
penguins |> 
    count(island) |>
    ggplot(aes(x=island, y=n)) +
    geom_col() +
    labs(title="a short title",
         subtitle=str_wrap("a really long title. kasjdhf ;isjdghf khg kajsxdhf khg alsidgf kjhg ljhags dfj hgkjahsdgfkjhg a  ljhsdgf ljhglsdjhfg", width=50))
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Instead use `ggtext` package - the element_textbox_simple will
automatically wrap text to fit whatever space is available.

``` r
library(ggtext)
```

``` r
penguins |> 
    count(island) |>
    ggplot(aes(x=island, y=n)) +
    geom_col() +
    labs(title="a short title",
         subtitle="a really long title. kasjdhf ;isjdghf khg kajsdhf khg alsidgf kjhg ljhags dfj hgkjahsdgfkjhg a  ljhsdgf ljhglsdjhfg") +
    theme(plot.subtitle = element_textbox_simple())
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

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

![](ggplot_demo_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

# discontinuous axes using `ggbreak` package

`ggbreak` package
[vignette](https://cran.r-project.org/web/packages/ggbreak/vignettes/ggbreak.html)

``` r
library(ggbreak) 
```

    ## ggbreak v0.1.2
    ## 
    ## If you use ggbreak in published research, please cite the following
    ## paper:
    ## 
    ## S Xu, M Chen, T Feng, L Zhan, L Zhou, G Yu. Use ggbreak to effectively
    ## utilize plotting space to deal with large datasets and outliers.
    ## Frontiers in Genetics. 2021, 12:774846. doi: 10.3389/fgene.2021.774846

``` r
library(patchwork)
```

there is a blank plot below, as well as the intended plots, but it
doesn’t appear when you knit to html or github_document

``` r
set.seed(2019-01-19)
d <- data.frame(x = 1:20,
   y = c(rnorm(5) + 4, rnorm(5) + 20, rnorm(5) + 5, rnorm(5) + 22)
)
 
p1 <- ggplot(d, aes(y, x)) + geom_col(orientation="y")
d2 <- data.frame(x = c(2, 18), y = c(7, 26), label = c("hello", "world"))
p2 <- p1 + scale_x_break(c(7, 17)) + 
  geom_text(aes(y, x, label=label), data=d2, hjust=1, colour = 'firebrick')  + 
  xlab(NULL) + ylab(NULL) + theme_minimal()

p1 + p2
```

![](ggplot_demo_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->
