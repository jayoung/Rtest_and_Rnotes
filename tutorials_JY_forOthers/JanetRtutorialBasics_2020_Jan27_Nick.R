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

# this is an example using a builtin dataset in R, that's often used as an example in help pages
myDat2 <- mtcars

# specify column(s):
myDat2[,"mpg"]
myDat2[,1:2]

# specify row(s):
myDat2[1,]
myDat2[20:30,]
myDat2[c(20,23,32),]

# add new columns:
myDat2[,"gearCarb"] <- myDat2[,"gear"] + myDat2[,"carb"] 
myDat2[,"gearCarbNorm"] <- (myDat2[,"gear"] + myDat2[,"carb"] ) / myDat2[,"qsec"] 

# write a table to file, to view later
write.table(myDat2, file="newTabDelimitedFile.txt", sep="\t")

##### means, max, min etc
mean(myDat2[,"mpg"])
median(myDat2[,"mpg"])
min(myDat2[,"mpg"])
max(myDat2[,"mpg"])
summary(myDat2[,"mpg"])
   ## if there are any NA datapoints (missing data) you need to do something like this:
max(myDat2[,"mpg"], na.rm=TRUE)


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

# a histogram
hist(myDat2[,"mpg"])



