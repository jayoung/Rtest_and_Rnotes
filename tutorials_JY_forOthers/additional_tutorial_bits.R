

save.image(file="expn.RData")

list.files() 
getwd()
setwd("Rdata")
setwd("..")
history(max.show=100)
write.table(arrayname,file="outfile.txt",sep="\t",row.names=FALSE,col.names=FALSE)

t(arrayname)
# t() transposes the matrix (swaps rows and columns)


# statistical tests

t.test()
wilcox.test()

chisq.test()

chisq.test()$obs
chisq.test()$exp


# strsplit, sapply
listresultingfromstrsplit <- strsplit(rawlanes[,"samplenorep"],"_")
allfirstelementsoflist <- sapply (listresultingfromstrsplit , "[" , 1)
allsecondelementsoflist <- sapply (listresultingfromstrsplit , "[" , 2)

# other apply functions

tapply()
lapply()



# random sampling
-----------
rbinom(100,1,0.35)
# samples 100 values from a binomial distribution of 1s and 0s, where 0.35 of the distribution are 1s.

rnorm()


# for loops
for (i in 1:10) {
    cat ("i ", i, "\n")
}

# if/else  conditionals

if (sum(col)==0) { 
    # do something 
    ngroups <- c(ngroups,"NA")
    next
} else {
    # do something else 
}

