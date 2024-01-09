## https://github.com/YuLab-SMU/treeio/issues/120

## I wondered whether there was a character/factor/numeric conversion error, but all these behave the same (wrong) way

## load package, get data, 
library(treeio)
data(bird.orders, package="ape")

#### add fake node labels - characters
b_nodeChars <- bird.orders
b_nodeChars$node.label <- paste("node",1:Nnode(b_nodeChars),sep="")
b_nodeChars_treedata <- as.treedata(b_nodeChars)

root(b_nodeChars, "Galliformes")
root(b_nodeChars_treedata, "Galliformes")


#### add fake node labels - integers
b_nodeIntegers <- bird.orders
b_nodeIntegers$node.label <- sample(100, Nnode(b_nodeIntegers))
b_nodeIntegers_treedata <- as.treedata(b_nodeIntegers)

root(b_nodeIntegers, "Galliformes")
root(b_nodeIntegers_treedata, "Galliformes")

#### add fake node labels - numeric
b_nodeNumeric <- bird.orders
b_nodeNumeric$node.label <-  sample(100, Nnode(b_nodeNumeric))/100
b_nodeNumeric_treedata <- as.treedata(b_nodeNumeric)

root(b_nodeNumeric, "Galliformes")
root(b_nodeNumeric_treedata, "Galliformes")




