require(graphics)
x  <- as.matrix(mtcars)
rc <- rainbow(nrow(x), start=0, end=.3)
cc <- rainbow(ncol(x), start=0, end=.3)
hv <- heatmap(x, col = cm.colors(256), scale="column",
              RowSideColors = rc, ColSideColors = cc, margin=c(5,10),
              xlab = "specification variables", ylab= "Car Models",
              main = "heatmap(<Mtcars data>, ..., scale = \"column\")")
str(hv) # the two re-ordering index vectors

## no column dendrogram (nor reordering) at all:
##heatmap(x, Colv = NA, col = cm.colors(256), scale="column",
  ##      RowSideColors = rc, margin=c(5,10),
    ##    xlab = "specification variables", ylab= "Car Models",
      ##  main = "heatmap(<Mtcars data>, ..., scale = \"column\")")

## "no nothing"
##heatmap(x, Rowv = NA, Colv = NA, scale="column",
  ##      main = "heatmap(*, NA, NA) ~= image(t(x))")

##round(Ca <- cor(attitude), 2)
##symnum(Ca) # simple graphic
##heatmap(Ca,             symm = TRUE, margin=c(6,6))# with reorder()
##heatmap(Ca, Rowv=FALSE, symm = TRUE, margin=c(6,6))# _NO_ reorder()

## For variable clustering, rather use distance based on cor():
##symnum( cU <- cor(USJudgeRatings) )

##hU <- heatmap(cU, Rowv = FALSE, symm = TRUE, col = topo.colors(16),
##             distfun = function(c) as.dist(1 - c), keep.dendro = TRUE)
## The Correlation matrix with same reordering:
##round(100 * cU[hU[[1]], hU[[2]]])
## The column dendrogram:
##str(hU$Colv)
