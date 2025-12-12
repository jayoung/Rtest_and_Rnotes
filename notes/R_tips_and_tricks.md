# R tips and tricks

To use a function that's not exported from a package (it's hidden), use three ::: symbols,  e.g. `seqLogo:::addLetter()`. If the function is exported then two :: are enough.


# Random bits and pieces

BEWARE: partial row name matching when using [ to subset.

Interpreting strings as variable names

One way:
```
x <- 42
eval(parse(text = "x"))
```

Another way:
```
x<-1:10
stringx <- "x"
get(stringx)
```



downloading files within R:
```
download.file(url, destfile="data-raw/name-of-file.xlsx")
```


Remove all variables whose names match the pattern "temp":
```
rm(list=ls(pat="temp"))
```

Rreading files from a .zip file:

```
## the first argument is the zip file, the second is the path that would result once the zip is unpacked
con <- unz("SRR10696710_fastqc.zip", "SRR10696710_fastqc/fastqc_data.txt")
dat <- scan(con, what="character")
```


`eapply` applies across all elements of an environment


str:  quick view on any object
dir.create
download.file

na.omit  (omits rows where values in any column are NA)
drop_na  (tidyverse version of na.omit)

Regular expressions. They're a bit different in R. The `glob2rx` function can help me convert a command-line style glob (e.g. "abc.*") into an R regexp (e.g. "^abc\\.")


### R does rounding weirdly!

See [here](https://psiaims.github.io/CAMIS/R/rounding.html)

"The `round()` function in Base R will round to the nearest whole number and ‘rounding to the even number’ when equidistant, meaning that exactly 12.5 rounds to the integer 12. Note that the janitor package in R contains a function `round_half_up()` that rounds away from zero. in this case it rounds to the nearest whole number and ‘away from zero’ or ‘rounding up’ when equidistant, meaning that exactly 12.5 rounds to the integer 13.""


### Insert an image into an Rmd document

(or any md document, I think)

A simple image:
```
![Caption for the picture.](/path/to/image.png)
```

Modify size (specify pixels, probably)
```
![Caption for the picture.](/path/to/image.png){#id .class width=30 height=20px}
```

Modify size by simple scaling
```
![Caption for the picture.](/path/to/image.png){#id .class width=50% height=50%}
```



To extract some list elements using their indices in a pipe, use `magrittr::extract()`. This is the equivalent of using `[`.   Example:

```
my_list |> magrittr::extract(3:5) 
```

The related function `magrittr::extract2()` is the equivalent of `[[`, so we would use it to extract a single element from a list. Example:

```
my_list |> magrittr::extract2(1) 
```

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



`unz()` function lets you read files within a zip file without even unpacking it. See also `gzfile()` and `bzfile()`


`nrow()` is a base R function. There is also a BiocGenerics function called `nrows()`, which I think is supposed to work on a list of matrix-like objects

If I have problems with package conflicts - different functions with the same name in different packages masking one another - I may be able to solve that using the `use()` function (available from R 4.4.0 onwards) to load functions from a library, rather than the `library()` function. . Using the `include.only` argument lets you load only the specified functions from a package rather than the whole package. Demo [here](https://erikgahner.dk/2025/use-use-in-r/).


`separate_longer_delim` is a nice function that lets us "uncollapse" rows when needed:
```
df <- data.frame(id = 1:3, x = c("x", "x y z", NA))
df
#   id     x
# 1  1     x
# 2  2 x y z
# 3  3  <NA>

df |> tidyr::separate_longer_delim(x, delim = " ")
#>   id    x
#> 1  1    x
#> 2  2    x
#> 3  2    y
#> 4  2    z
#> 5  3 <NA>
```


# getting help

see
http://journal.r-project.org/archive/2009-2/RJournal_2009-2_Graves~et~al.pdf

```
help(plot)
help.search("plot")

help(package="BiocInstaller")
will list the functions

library(sos) 
PL <- findFn("Petal.Length")  
#### seems like it goes through lots of R docs, finds help pages for any function that mentions Petal.Length. When I then call PL, it shoots up an html page with links to all of those

RSiteSearch("Petal.Length")  
### shoots up a web page with links to lots of things mentioning Petal.Length. From utils package. 

library(Biobase)
openVignette(package="seqLogo") #### openVignette is a Biobase function

vignette()  
   lists all vignettes
vignette(package="grid")
   lists grid's vignettes
vignette("saveload")  
   to load a specific vignette (this one is from the grid package)
vignette("rotated", package="grid")
   might need to do this to load vignettes from specific packages??

or: 
openVignette("ShortRead")

library(help = grid)
   lists all functions in a package and vignettes
   
v1 <- vignette("grid")
print(v1)  ### brings up the pdf
edit (v1)  ### brings up the associated R code

v1 <- vignette("KEGGgraph")
print(v1)  ### brings up the pdf
edit (v1)  ### brings up the associated R code
```

## Debugging


Three options:    
- `browser()` (place inside a function, temporarily)    
- `debug(myFunction)` plus `undebug(myFunction)` (the Rstudio console window has a 'stop' button to exit debugging)   
- `debugonce()`    
See (`explore_debugging_functions.R`)[Rscripts/explore_debugging_functions.R] for details.

[A tutorial](https://data-flair.training/blogs/debugging-in-r-programming/) on debigging

## Pipes:

I was using `magrittr`'s pipe in my code:   `%>%`.   

I recently switched MOST code in some of my repos to the native R pipe: `|>` (introduced in R 4.1.0).

Can insert a pipe in Rstudio using Command-shift-m on the mac. There is a setting in R studio to tell it which of those two styles of pipe you want to use.

There are situations where the native pipe won't work and we need to use  `%>%`. More info [here](https://tidyverse.org/blog/2023/04/base-vs-magrittr-pipe/). Examples: 
- when we use the `.` to represent the incoming data as an argument to one of the functions.  
- if we don't include the parentheses after the function name, e.g. `my_dat |> nrow()` works but `my_dat |> nrow` doesn't. (and `my_dat %>% nrow` does work)

`nrow()` has caused me trouble in pipes. A better alternative is the `tibble::rowid_to_column()` function, e.g. `iris |> rowid_to_column(var="my_row_index") |> head()`

There are other [more complicated pipe types](https://cran.r-project.org/web/packages/magrittr/vignettes/magrittr.html), including this pipe: `%<>%`. But I doubt I'll use them because it would make the code quite unreadable. 


# reprex

To create reproducible code + output we can use `reprex`

First, we write the code the demonstrates the problem

Load the library:
```
library(reprex)
```

Copy the code you want to the clipboard, and enter `reprex()` in the R Console. In RStudio, you’ll see a preview of your rendered reprex, but it is also now on the clipboard ready to paste.

Or, write the code, and then put it in a reprex block:
```
reprex::reprex({
    x <- 1:4
    y <- 2:5
    x + y
})
```



# Memory usage 

```
gc - looks useful
object.size
mem.limits
?Memory
Rprof
summaryRprof
``` 


Rprofmem and tracemem (must be enabled at compile time)

memory.profile - doesn't look useful


# use secondary y-axis in base R

http://www.r-bloggers.com/r-graph-with-two-y-axes/

```
x <- 1:5
y1 <- rnorm(5)
y2 <- rnorm(5,20)

plot(x,y1,type="l",col="red")
par(new=TRUE)
plot(x, y2,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("y2",side=4,line=3)
legend("topleft",col=c("red","blue"),lty=1,legend=c("y1","y2"))
```


# Rstudio 


## Rstudio tricks

"Reindenting your code only shifts things around horizontally. If you want more powerful code reformatting, try using “Code > Reformat Code” (or use ⌘⇧A on macOS or ctrl + shift + A on Windows). It’s a more aggressive form of reformatting that will add extra line breaks and other things to make the code more readable."

In settings, there's an option to turn on rainbow parentheses to help see pairings.


## Rstudio keyboard shortcuts

highlight a function, press function-F1, and it brings up the help page

shift-command-M types a pipe (there's a setting for whether you want that to be `%>%` or `|>`)

tab adds an indent to one or more selected lines of code, shift+tab removes one indent

(see tips [here](https://rfortherestofus.com/2023/11/rstudio-hotkeys))

## Rstudio 'snippets'

[Snippets](https://rfortherestofus.com/2024/04/snippets-rstudio):

- e.g. Example: type `fun` and press the `tab` key, and R provides the skeleton of a new function
- press tab again and you move within the snippet to the next piece you might fill in. Shift-tab also does something (not sure exactly what)
- `Tools menu - Edit code snippets` shows what snippets are available
- Markdown snippets are also useful (e.g. place an image). Here we need to do `shift-tab` to activate.  E.g. `r-shift-tab` inserts an R code chunk, if we do it from within a markdown area (i.e. outside an existing R code chunk)

To see all snippets:  Tools - Edit Code Snippets


## Rstudio speed issues

sometimes have trouble with Rstudio acting super slow and weird. Seems to be with things stored on network drive, particularly when using Rprojects in dirs I'm syncing to github.

Advice:
https://community.rstudio.com/t/rstudio-startup-time-40-seconds/26174/8?u=kevinushey
disable a bunch of stuff in Rstudio preferences:
- disable git/svn version control
- disable everything under code-diagnostics



# Positron

An IDE for R or python.  Perhaps a future replacement for Rstudio but as of 2025 it's not able to interact with data stored on the Hutch servers, although I could use it on my Mac, I think.

Its' data viewer looks nice (`View()`) - includes little graphs


I am trying it on the work laptop (Dec 2025). Some notes from the walkthrough it gives on migrating from Rstudio:

- "Positron doesn't have an exact equivalent to RStudio Projects, but the concept in Positron that is most analogous to an RStudio Project is a workspace. You can read more in the VS Code documentation about what exactly a workspace is, but in general think of a workspace as about the same thing as a folder, which is about the same thing as an RStudio Project."  More info on [how to think about Rprojects in Positron](https://positron.posit.co/migrate-rstudio-rproj.html)
- "Air" is a formatting tool for R code
- "Databot" is a Positron extension that allows you to use AI, e.g. via the Claude LLM, to do data analysis, including saving the R code for it. [My notes](https://github.com/jayoung/MalikLab_bioinformaticsResources/blob/main/janets_NOTES_forMyself/programming_and_statistics/AI_notes.md#databot-in-positron-for-r-based-data-analysis) (in a different repo).


# showing colors 

```
totalColors <- length(colors())
# 657 of them

plotCols <- function(firstColorIndex, lastColorIndex) {
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),  
         axes = FALSE, xlab = "", ylab = "")
    
    numRows <- 42
    numColumns <- 8
    textXoffset <- 0.1 * 1/numColumns
    textYoffset <- 0.375 * 1/numRows
    myColors <- colors()[seq(firstColorIndex,lastColorIndex)]
    
    rect( rep((0:(numColumns - 1)/numColumns),numRows) ,
          sort(rep((0:(numRows - 1)/numRows),numColumns),decreasing=T),
          rep((1:numColumns/numColumns),numRows) , 
          sort(rep((1:numRows/numRows),numColumns),decreasing=T), 
              border = "light gray", 
              col=myColors)
    
    text( rep((0:(numColumns - 1)/numColumns),numRows) + textXoffset , 
          sort(rep((0:(numRows - 1)/numRows),numColumns),decreasing=T) + textYoffset ,  
          myColors, adj=0,
        cex=0.5)

}

pdf(file="RcolorDemo.pdf", width=7,height=11)
par(mar=c(0,0,0,0))
plotCols(1, 336)
plotCols(337, totalColors)
dev.off()
```



# Using specific fonts

older notes on how to use Arial font, in R: (modified from [PLoS Genet instructions](http://www.plosgenetics.org/static/figureGuidelines.action#arialR)

```
postscript(file="try.ps", horizontal=F,onefile=F,width=4, height=4,
   family=c("/home/jayoung/fontFiles/arial.afm",
            "/home/jayoung/fontFiles/arialbd.afm",
            "/home/jayoung/fontFiles/ariali.afm",
            "/home/jayoung/fontFiles/arialbi.afm"),  pointsize=12)
hist(rnorm(100))
dev.off()
```
For PLoS submissions, it seems I have to save from R as a postscript, to make sure fonts aren't screwed up (when I save as pdf, some points in graphs use an unidentified font, which Illustrator converts to AdobePiStd, which PLoS website doesn't recongize). Then I use AI to save as eps (make sure I choose encapsulated fonts)

More advice on fonts [here](http://r.789695.n4.nabble.com/pdf-device-uses-fonts-to-represent-points-data-alteration-td836140.html)

To make the figures look the right size for my PLoS Genet submission, I also changed three lines in the eps file using a text editor (if I don't do that, they re-scale to fill the whole page): 
```
%%BoundingBox: 0 0 612 792
%%HiResBoundingBox: 0 0 612 792
%%CropBox: 0 0 612 792
```




## R infrastructure

Hidden files that are read when R starts up:

```
.libPaths()
~/.Rprofile
./.Rproj.user
./*Rproj
./.Rhistory
.Renviron
```

You can [run R without any startup files](https://rstats.wtf/r-startup.html) by using the --vanilla argument when starting R


`.Rprofile` gets read when I start up R, so I can put functions/options I use regularly in there.  If there's one in the project directory (local), that one will be read. If not, it'll look for one in `~` (global).  If you want to load BOTH local and global `.Rprofile` files, you need to include a line in the local file that looks something like this: `source("/home/jayoung/.Rprofile")`.  


To investigate environmental variables from within R:
```
Sys.getenv("LD_LIBRARY_PATH")
temp <- Sys.getenv()
names(temp)
temp[grep ("btrask",temp)]
temp[grep ("lib_linux",temp)]

Sys.getenv("R_RD4PDF")
```
Note - paths will probably be different within R than outside it.

There is also a file called `Renviron` installed with R that will generate some of these environmental variables.  Example on gizmo/rhino: `/app/software/R/4.4.1-gfbf-2023b/lib/R/etc/Renviron`



On a Mac, R packages go here - `/Library/Frameworks/R.framework/Versions` - in subdirectories by version. After installing new R, can delete old packages to save disk space