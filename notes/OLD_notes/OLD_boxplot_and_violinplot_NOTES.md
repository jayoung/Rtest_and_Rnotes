

# boxplots:

Base R:
```
boxplot(omega~nonconscode,
       data=data98yesno,
       xlab="number of non-consensus residues per pair", 
       ylab = "omega", main="noncons98, LeoMR05nonpseud, omega")
```
omega and nonconscode are column titles in the data matrix - this will plot the omegas classified by the tag in nonconscode
xlab is x axis label, main is title

let's say I want to make a boxplot of a dataframe (called df), where each boxplot item is a column of the data.frame:

```
boxplot(list(x=x,y=y))

boxplot(df, (rep(colnames(df),each=dim(df)[1])) )
```

# Test data

Make some test data: a bimodal distribution (x1 or x2), with means at 0 (2/3 of the datapoints) and 5 (1/3 of the datapoints).  Also, make a factor (f1 or f2) that divides the data into 5 groups. 

```
mySmallNumber <- 100  ### my dataset will be 1.5 times this big:
x1 <- c(rnorm(mySmallNumber), rnorm(mySmallNumber/2,5))
f1 <- factor(rep(1:5, (1.5*mySmallNumber)/5))
df1 <- data.frame(x1=x1, f1=f1)

myBigNumber <- 100000  ### my dataset will be 1.5 times this big:
x2 <- c(rnorm(myBigNumber), rnorm(myBigNumber/2,5))
f2 <- factor(rep(1:5, (1.5*myBigNumber)/5))
df2 <- data.frame(x2=x2, f2=f2)
```

# Violin plots

I was investigating how long some of these functions took.  Some were very slow.

```
# simply plot all the data as one group:
system.time(PlotViolin(x2))
   ## 90 seconds!!
system.time(vioplot(x2))
    ## < 1 second but I don't like the way it chooses bandwidth - too coarse-grained
system.time(PlotViolin(x2, bw="nrd0"))
   ## much quicker!  and, I like the bandwidth choice. 
```

ggplot2
```
ggplot2.violinplot(data=x2)

ggplot2.violinplot(data=x2, addMean=TRUE, meanPointShape=23, meanPointSize=3, meanPointColor="black", meanPointFill="blue" )

violinplot(x2,col="brown")
```

Subsetting the data: I previous had not figured out with any of these functions how to zoom in on a portion of the distrubution - I can subset the data that goes into the distribution plot, but I cannot figure out how to just truncate the portion that is displayed.  Now I know it's `+coord_cartesian()`


```
ggplot2.violinplot(data=x2, ylim=c(0,5))
# Warning message:
# Removed 75079 rows containing non-finite values (stat_ydensity). 

violinplot(x2~f2,col="brown", ylim=c(0,5) )

#### don't know how to do it with PlotViolin - neither of these works:
#  PlotViolin(x2, bw="nrd0", subset=x2<5 ) 
#  PlotViolin(x2 ~ f2, bw="nrd0", xlab="x axis", ylab="y axis", ylim=c(0,5) ) 
```

plot the data by groups: 
```

system.time(PlotViolin(x2 ~ f2, bw="nrd0"))

ggplot2.violinplot(data=df2, xName="f2",yName="x2",addMean=TRUE)

violinplot(x2~f2,col="brown")
```

For vioplot, need a special hacked version of the function in order to pass in a list as input: got this from http://stackoverflow.com/questions/22410606/violin-plot-with-list-input:

```
vioplot2( split(x2, f2) )

vioplot2<-function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
    horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
    lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
    at, add = FALSE, wex = 1, drawRect = TRUE) 
{
    if(!is.list(x)){
        datas <- list(x, ...)
    } else{
        datas<-x
    }
    n <- length(datas)
    if (missing(at)) 
        at <- 1:n
    upper <- vector(mode = "numeric", length = n)
    lower <- vector(mode = "numeric", length = n)
    q1 <- vector(mode = "numeric", length = n)
    q3 <- vector(mode = "numeric", length = n)
    med <- vector(mode = "numeric", length = n)
    base <- vector(mode = "list", length = n)
    height <- vector(mode = "list", length = n)
    baserange <- c(Inf, -Inf)
    args <- list(display = "none")
    if (!(is.null(h))) 
        args <- c(args, h = h)
    for (i in 1:n) {
        data <- datas[[i]]
        data.min <- min(data)
        data.max <- max(data)
        q1[i] <- quantile(data, 0.25)
        q3[i] <- quantile(data, 0.75)
        med[i] <- median(data)
        iqd <- q3[i] - q1[i]
        upper[i] <- min(q3[i] + range * iqd, data.max)
        lower[i] <- max(q1[i] - range * iqd, data.min)
        est.xlim <- c(min(lower[i], data.min), max(upper[i], 
            data.max))
        smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
            args))
        hscale <- 0.4/max(smout$estimate) * wex
        base[[i]] <- smout$eval.points
        height[[i]] <- smout$estimate * hscale
        t <- range(base[[i]])
        baserange[1] <- min(baserange[1], t[1])
        baserange[2] <- max(baserange[2], t[2])
    }
    if (!add) {
        xlim <- if (n == 1) 
            at + c(-0.5, 0.5)
        else range(at) + min(diff(at))/2 * c(-1, 1)
        if (is.null(ylim)) {
            ylim <- baserange
        }
    }
    if (is.null(names)) {
        label <- 1:n
    }
    else {
        label <- names
    }
    boxwidth <- 0.05 * wex
    if (!add) 
        plot.new()
    if (!horizontal) {
        if (!add) {
            plot.window(xlim = xlim, ylim = ylim)
            axis(2)
            axis(1, at = at, label = label)
        }
        box()
        for (i in 1:n) {
            polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
                c(base[[i]], rev(base[[i]])), col = col, border = border, 
                lty = lty, lwd = lwd)
            if (drawRect) {
                lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
                  lty = lty)
                rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
                  q3[i], col = rectCol)
                points(at[i], med[i], pch = pchMed, col = colMed)
            }
        }
    }
    else {
        if (!add) {
            plot.window(xlim = ylim, ylim = xlim)
            axis(1)
            axis(2, at = at, label = label)
        }
        box()
        for (i in 1:n) {
            polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                rev(at[i] + height[[i]])), col = col, border = border, 
                lty = lty, lwd = lwd)
            if (drawRect) {
                lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
                  lty = lty)
                rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
                  boxwidth/2, col = rectCol)
                points(med[i], at[i], pch = pchMed, col = colMed)
            }
        }
    }
    invisible(list(upper = upper, lower = lower, median = med, 
        q1 = q1, q3 = q3))
}
```

