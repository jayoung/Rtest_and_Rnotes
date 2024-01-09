# Rtest_and_Rnotes
my R playground and notes

rhino location: `~/FH_fast_storage/git_more_repos/Rtest_and_Rnotes`


## Current learning
[Tidyverse style guide](https://style.tidyverse.org/syntax.html).  I got up to the end of the 'syntax' section.  Perhaps see also the [Advanced R style guide](http://adv-r.had.co.nz/Style.html).

## Resources

Nick Tierney's (mostly) [rstats blog](https://www.njtierney.com)



### haven't looked yet: future

[The Elements of Data Analytic Style](https://leanpub.com/datastyle) by Jeff Leek (a Leanpub book)


## Things I've learned

### Debugging

Three options:    
- `browser()` (place inside a function, temporarily)    
- `debug(myFunction)` plus `undebug(myFunction)`    
- `debugonce()`    
See (`explore_debugging_functions.R`)[Rscripts/explore_debugging_functions.R] for details.

### Misc 

The 'embracing' operator (`{{ }}`), and unquoting using !! and !!! - see [`testCode.R`](Rscripts/testCode.R) for details.

& versus && (and | versus ||):  use the short form for bitwise operation on vectors. Use the long form when we want a single TRUE/FALSE answer.   `any` and `all` functions run OR and AND on all elements of a vector

```
x <- c(TRUE,TRUE,FALSE,FALSE)
y <- c(TRUE,FALSE,TRUE,FALSE)
x & y
# x && y # this is no good!
any(x)  ## TRUE
all(x)  ## FALSE
```

The `switch` function - a multiway `if` statement, I think?
```
centre <- function(x, type) {
  switch(type,
         mean = mean(x),
         median = median(x),
         trimmed = mean(x, trim = .1))
}
x <- rcauchy(10)
centre(x, "mean")
centre(x, "median")
centre(x, "trimmed")
```

In `switch`, if there's an empty argument, it 'falls-through' to the next thing (e.g. here, `myFunc("a")` returns the same thing as `myFunc("b")`).
Not also that we can add `call. = FALSE` to a `stop` statement to modify the error message that'll be produced
```
myFunc <- function(x) {
  switch(x, 
    a = ,
    b = 1, 
    c = 2,
    stop("Unknown `x`", call. = FALSE)
  )
}
```


## Useful packages

### phylogenetics

ape

[ggtree](https://yulab-smu.top/treedata-book/)  (also tidytree and treeio)

[ggtreeExtra](https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html).  ggtree can use geom_facet to align associated graphs to the tree but it only works with rectangular, roundrect, ellipse and slanted layouts. ggtreeExtra allows graphs on a tree in rectangular, circular, fan and radial layouts

in ~/domesticated_capsid/Rreports/RTL3_frameshift_plots_v2_aln28.Rmd I got a tree of >5000 mammal species from Upham publication, and extracted the species I want