########## defining variables

x1 <- 3
x2 <- 10
x3 <- c(3,10)
x4 <- 1:10
x5 <- "apple"
x6 <- c("apple", "pear", "orange")
x7 <- c("apple", "apple", "apple", "pear", "orange")

## tips: don't use variable names that might be function names in R (e.g. table, mean) - try something that's more likely to be unique (e.g. myTable, myMean)

########## a few useful functions

length(x6)
mean(x4)
median(x4)
summary(x4)
table(x7)
range(x4)

# subsetting:
x7[2]
x7[2:3]
x7[c(1,4,5)]

########## object classes.  each variable/object has a class, and R is picky about what it can do with each class.

class(x1)
class(x7)

## typical R object classes:
## integer, numeric, factor, character, data.frame, matrix, list

########## Rstudio - I haven't used this but I think many people find it helpful. (e.g. Tera?)

########## getting help

?mean
help("mean")
help.search("average")
     # the examples at the bottom of each help page are very useful

# google search - include r-help, e.g. :  search for    r-help median

# https://cran.r-project.org/
#    manuals - an introduction to R
#    online tutorials???

########## working with tables - very common usage.  I find data.frame objects more intuitive than matrix objects (more flexible)

?read.delim
myDat <- read.delim("myTabDelimitedFile.txt")
myDat <- read.delim("mySpaceDelimitedFile.txt", header=FALSE, sep=" ")

# this was an example using one of Michelle's files
myDat <- read.delim("Permissivepool1_chrI.txt", header=FALSE, sep=" ") 

head(myDat)
tail(myDat)
dim(myDat)
colnames(myDat) <- c("pos","allele1","allele2")

# specify column(s):
myDat[,"pos"]
myDat[,1:2]

# specify row(s):
myDat[1,]
myDat[20:30,]
myDat[c(20,23,32),]

# add a new column:
myDat[,"totalReads"] <- myDat[,"allele1"] + myDat[,"allele2"] 
myDat[,"allele1freq"] <- myDat[,"allele1"] / (myDat[,"allele1"] + myDat[,"allele2"])

write.table(myDat, file="newTabDelimitedFile.txt", sep="\t")


########## making plots

plot (1:10, 1:10)
plot (1:10, 1:10, pch=19, cex=0.5, xlab="my x axis", ylab="my y axis", main="my plot title")
plot (1:10, 1:10, pch=19, cex=0.5, xlab="my x axis", ylab="my y axis", main="my plot title", xlim=c(0,5))

?plot
?par

# several plots per page:
X11(height=11, width=7)
par(mfrow=c(2,1))
plot (1:10, 1:10)
plot (1:10, 1:10, pch=19, xlab="my x axis", ylab="my y axis")

plot(density(myDat[,"totalReads"]))

plot(density(myDat[,"totalReads"], from=0,to=200))

plot(density(myDat[,"totalReads"], from=0,to=200))
lines(density(myDatPermissive[,"totalReads"], from=0,to=200), col="red")


########## splitting a data frame into a list of small dataframes

table(table(myDat[,"cloneSequence"]))

myDatSplit <- split(myDat, myDat[,"cloneSequence"])

## subsetting for lists:
myDatSplit[[1]]  ## double-brackets for a single element
myDatSplit[1:3]  ## single brackets for a few elements (I don't know why)

class(myDatSplit)
class(myDatSplit[[1]])

##### applying a function to all the things in the list:  

# lapply returns a list
lapply(myDatSplit, dim)

# sapply shows the results in a tidier way, if possible
sapply(myDatSplit, dim)

# sometimes the function is more complicated and we need to wrap it up in parentheses, or even define it as a named function that we can reuse:

sapply( myDatSplit, function(x) {
    answer <- mean(x[,"measurement"])
    return(answer)
} )


myNamedFunction <- function(x) {
    answer <- mean(x[,"measurement"])
    return(answer)
} )
sapply(myDatSplit,myNamedFunction)


## maybe we choose and keep only the row with the highest measurement, and we add a column that states how many rows there were to begin with
myDatSplitBest <- lapply(myDatSplit, function(x) {
    y <- x[ which(x[,"measurement"] == max(x[,"measurement"]))[1], ]
    y[,"numClones"] <- dim(x)[1]
    return(y)
})

#### re-combining a list of data.frames
myDatSplitRejoined <- do.call("cbind",myDatSplitBest) 


##### for simple uses, tapply does the splitting AND the applying:
tapply( myDat[,"measurement"], myDat[,"cloneSequence"], mean)  # get mean according for each clone sequence
tapply( myDat[,"measurement"], myDat[,"cloneSequence"], sum)

########## Bioconductor

# there are various packages of functions to work with typical biological sequence data. The Bioconductor project hosts many packages, and there are also regular R packages (not specifically for biological data)

# to load any package:
library(VariantAnnotation)

# a few Bioconductor object classes:
# GRanges, CollapsedVCF, Rle

# many Bioconductor packages have "vignettes" demonstrating how to use the package:
openVignette("VariantAnnotation")
(or available online at bioconductor website)

