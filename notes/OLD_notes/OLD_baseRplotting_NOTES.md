# Big picture

All these notes are really old. Now I would always use ggplot2 - can't imagine going back.

## XY scatterplots:

```
plot(alldS,alldN,type="n",xlim=c(0,0.5),ylim=c(0,0.2),xlab="dS",ylab="dN",main="dN dS, all non-pseud pairs, by Zhang MOE expression status")
points(MOEenriched[,"dS"],MOEenriched[,"dN"],pch=20,cex=.5,col="red")
points(MOEenriched_marginal[,"dS"],MOEenriched_marginal[,"dN"],pch=20,cex=.5,col="green")
points(MOEnotenriched[,"dS"],MOEnotenriched[,"dN"],pch=20,cex=.5,col="black")

leg.txt<-c("MOEenriched","MOEenrichedmarginal","MOEnotenriched")
legend("topright",leg.txt,pch=20,col=c("red","green","black"))
```

## histograms

```
hist(recS,xlab="RecombSHRSPxBN",br=20,main="Recombination rate all intervals")
```
recS is an array
br=20 gives 20 categories
main is title
xlab is x label

to plot two histograms side by side, the data must be in list form, and use multhist from plotrix package:
```
require(plotrix)
mylist<-list(OR=OR[,2],V1R=V1R[,2])
cols<-grey.colors(2)
multhist(mylist,col=cols,main="Genomic clusters, using 0.5Mb cutoff, num genes")
leg.txt<-c("OR","V1R")
legend("topright",leg.txt,col=cols,fill=cols)
```

## 2D density plot

```
library(KernSmooth)
x<-cbind(ORs$Testis_1,ORs$OE_M0)
est<-bkde2D(x,bandwidth=c(10,10))
contour(est$x1,est$x2,est$fhat,xlim=c(-500,500),ylim=c(-500,500))
```
(makes contour plot of density function of points in the two columns of ORs data table labelled Testis_1 and OE_M0. change the bandwidth to alter smoothing)

or, using smoothScatter:
```
pairs(exprs[,allsamples], panel=function(...) {par(new=TRUE);smoothScatter(...,nrpoints=0)})
```

or, without smoothing (but this is inefficient with large data)
```
Lab.palette <- colorRampPalette(c("blue", "orange", "red"), space = "Lab")
dcols <- densCols(x,y,colramp = Lab.palette)
graphics::plot(x, y, col = dcols, pch=20, main = "n = 1000")
```

## mosaic plots for categorical data

```
plot(familycounts,
     xlab="Family",ylab="OrthoFuller",
     main="OrthoFuller by Family",
     color=rainbow(8),las=1)
```
(familycounts is a table of counts, las=1 just specifies category labels always horizontal, color allows you to specify the palette used - this case there were 8 categories, so palette includes 8 colors)

## pairs plots

```
upPan<- function (...) {
    points(...,col="darkblue")
    abline(a=0,b=1,col="red")
}
lowPan<-function(x,y,...) {
    text(# determine co-ordinates for the midpoint of the plot
         mean(par("usr")[1:2]), mean(par("usr")[3:4]), 
         signif(cor(x,y,use="complete.obs"),2),
         cex=2)
}
pairs(meths2[,braincols],pch=".", upper.panel=upPan, lower.panel=lowPan)
```

## Other

`venn()` for >3 way venn diags


Once in a while I end up with a REALLY big pdf file that I want to convert to png perhaps to make it more manageable size to print. I can do this from the linux command line:
```
convert -density 600x600 a.pdf a.png 
```
