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

########## object classes

class(x1)
class(x7)

typical R object classes:
integer, numeric, factor, character, data.frame, matrix, list

########## Bioconductor

there are various packages of functions to work with typical biological sequence data. The Bioconductor project hosts many packages, and there are also regular R packages (not specifically for biological data)

to load any package:
library(VariantAnnotation)

a few Bioconductor object classes:
GRanges, CollapsedVCF, Rle

########## Rstudio - I haven't used this but I think many people find it helpful. Tera can show you (?)

########## getting help

?mean
help("mean")
help.search("average")
     look at the examples at the bottom of each help page

google search:  "r-help mean"

many Bioconductor packages have "vignettes" demonstrating how to use the package:
openVignette("VariantAnnotation")


https://cran.r-project.org/
    manuals - an introduction to R
    online tutorials???

########## working with tables

?read.delim
myDat <- read.delim("myTabDelimitedFile.txt")
myDat <- read.delim("myTabDelimitedFile.txt", header=FALSE, sep=" ")

myDat <- read.delim("Permissivepool1_chrI.txt", header=FALSE, sep=" ") 

head(myDat)
tail(myDat)
dim(myDat)
colnames(myDat) <- c("pos","allele1","allele2")

specify column(x):
myDat[,"pos"]
myDat[,1:2]

specify row(s):
myDat[1,]
myDat[20:30,]
myDat[c(20,23,32),]

add a new column:
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


########## loops, list objects

list.files()

myFiles <- list.files(pattern="Permissivepool")

# simple for loop, to populate a list: 
datEachChr <- list()
for (thisFile in myFiles) {
    datEachChr[[thisFile]] <- read.delim(thisFile, header=FALSE, sep=" ") 
}
class(datEachChr)
class(datEachChr[[1]])

head(datEachChr[[1]])

names(datEachChr)
names(datEachChr) <- myFiles

datEachChr[["Permissivepool1_chrXVI.txt"]]

# a more elegant way:  use lapply (list apply)
datEachChr <- lapply( myFiles, function(x) {
    y <- read.delim(x, header=FALSE, sep=" ") 
    return(y)
})

# more complicated version of lapply
datEachChr <- lapply( myFiles, function(x) {
    y <- read.delim(x, header=FALSE, sep=" ") 
    colnames(y) <- c("pos","allele1","allele2")
    y[,"totalReads"] <- y[,"allele1"] + y[,"allele2"] 
    y[,"allele1freq"] <- y[,"allele1"] / (y[,"allele1"] + y[,"allele2"])
    return(y)
})
names(datEachChr) <- gsub("Permissivepool1_","",myFiles)
names(datEachChr) <- gsub(".txt","",names(datEachChr) )


lapply( datEachChr, head)


# also read in restrictive pool data
myFiles2 <- list.files(pattern="Restrictivepool")
datEachChrRestrictive <- lapply( myFiles2, function(x) {
    y <- read.delim(x, header=FALSE, sep=" ") 
    colnames(y) <- c("pos","allele1","allele2")
    y[,"totalReads"] <- y[,"allele1"] + y[,"allele2"] 
    y[,"allele1freq"] <- y[,"allele1"] / (y[,"allele1"] + y[,"allele2"])
    return(y)
})
names(datEachChrRestrictive) <- gsub("Restrictivepool1_","",myFiles2)
names(datEachChrRestrictive) <- gsub(".txt","",names(datEachChrRestrictive) )

# plot allele freqs for each chr:
pdf(file="myPlots.pdf", width=7, height=11)
par(mfrow=c(6,1))
for (thisChr in names(datEachChr) ) {
    thisDat <- datEachChr[[thisChr]]
    thisDat2 <- datEachChrRestrictive[[thisChr]]
    plot( thisDat[,"pos"]/1000, thisDat[,"allele1freq"], pch=19, cex=0.25, xlab="position (kb)", ylab="allele1freq", main=thisChr)
    points( thisDat2[,"pos"]/1000, thisDat2[,"allele1freq"], pch=19, cex=0.25, col="red")
    abline(h=0.5, col="gray", lty=2)
}
dev.off()

# plot coverage for each chr:
pdf(file="myPlots_coverage.pdf", width=7, height=11)
par(mfrow=c(6,1))
for (thisChr in names(datEachChr) ) {
    thisDat <- datEachChr[[thisChr]]
    thisDat2 <- datEachChrRestrictive[[thisChr]]
    plot( thisDat[,"pos"]/1000, thisDat[,"totalReads"], pch=19, cex=0.25, xlab="position (kb)", ylab="totalReads", main=thisChr, ylim=c(0,300))
    points( thisDat2[,"pos"]/1000, thisDat2[,"totalReads"], pch=19, cex=0.25, col="red")
}
dev.off()


# a way to combine a list of data.frames:
allPermissive <- do.call( "rbind", datEachChr)
allRestrictive <- do.call( "rbind", datEachChrRestrictive)

# plot coverage for all SNPs, each pool
plot(density(allRestrictive[,"totalReads"], from=0, to=200), xlab="total reads / SNP", ylab="relative frequency", main="coverage per SNP. black=restrictive pool, red=permissive pool" )
lines (density(allPermissive[,"totalReads"], from=0, to=200), col="red" )

# plot allele freq for all SNPs, each pool
plot(density(allPermissive[,"allele1freq"]), col="red" , xlab="allele1freq", ylab="relative frequency", main="allele1freq. black=restrictive pool, red=permissive pool" )
lines (density(allRestrictive[,"allele1freq"]))
abline(v=0.5, col="gray", lty=2)

# plot coverage versus allele freq for all SNPs, each pool
X11(height=11,width=7)
par(mfrow=c(2,1))

smoothScatter(  allPermissive[,"totalReads"], allPermissive[,"allele1freq"], xlab="coverage", ylab="allele1freq", pch=19, cex=0.25, main="permissive pool" )
abline(v=0.5, col="gray", lty=2)

smoothScatter(  allRestrictive[,"totalReads"], allRestrictive[,"allele1freq"], xlab="coverage", ylab="allele1freq", pch=19, cex=0.25, main="restrictive pool" )
abline(v=0.5, col="gray", lty=2)

