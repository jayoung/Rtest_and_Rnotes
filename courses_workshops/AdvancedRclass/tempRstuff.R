#### install this at some point: http://www.tug.org/mactex/
#### or maybe not because it takes a lot of disk space

source("http://bioconductor.org/course-packages/courseInstall.R")
courseInstall("AdvancedR2011")

library("AdvancedR2011Data")
library("AdvancedR2011")

browseVignettes("AdvancedR2011")
browseVignettes("AdvancedR2011Data")

### http://bioconductor.org/course-packages/snapshots/StudentGWAS_0.4.0.tar.gz

######### session 1 code

reduce <- function(lst) {
    n <- sum(sapply(lst, "[[", "N"))
    het <- Reduce("+",lapply(lst,"[[","Het"))
    het / n
}

get <- function(con) {
    f2(con, nrows=100,ncols=.ncols, keep=1:100)
}

heterozygosity <- function (chunk, ...) {
    list(N=nrow(chunk),Het=colSums(chunk==2))
}


######### session 2: SQL

######### get file names of the csv files containing the data, and read in that data
list.files( system.file("extdata",package="AdvancedR2011Data") )
sub <- list.files( system.file("extdata",package="AdvancedR2011Data"), full.names=T )[10]
snps <- list.files( system.file("extdata",package="AdvancedR2011Data"), full.names=T )[7]
sub_dat <- read.csv(sub)
snps_dat <- read.csv(snps)

###### make a database
library(RSQLite)
drv <- SQLite()
con <- dbConnect(drv, dbname="metadata.sqlite")

#### make two tables in that database
dbGetQuery(con, "CREATE Table subjects (id INTEGER, subject_id TEXT, case_control TEXT)")
dbGetQuery(con, "CREATE Table snps (id INTEGER, seqnames TEXT, start INTEGER, gene_id TEXT, snp_id TEXT)")

### populate the tables
sql <- "INSERT INTO subjects VALUES ($id, $subject_id, $case_control)"
dbBeginTransaction(con)
dbGetPreparedQuery(con, sql, bind.data = sub_dat)
dbCommit(con)

sql <- "INSERT INTO snps VALUES ($id, $seqnames, $start, $gene_id, $snp_id)"
dbBeginTransaction(con)
dbGetPreparedQuery(con, sql, bind.data = snps_dat)
dbCommit(con)

##### check on database contents
dbListTables(con)
dbListFields(con,"subjects")
dbListFields(con,"snps")

#### create indices for those tables
sql <- "CREATE INDEX subId ON subjects(id)"
dbGetQuery(con, sql)
sql <- "CREATE INDEX snpId ON snps(id)"
dbGetQuery(con, sql)

### finish up
dbDisconnect(con)


### reopen that database
library(RSQLite)
drv <- SQLite()
mycon <- dbConnect(drv, dbname="~/temp/StudentGWAS/inst/extdata/metadata.sqlite")

getSnps <- function (dbcon) {
    sql <- "SELECT * FROM snps"
    return(dbGetQuery(dbcon, sql))
}
head(getSnps(mycon))
class(getSnps(mycon))
dim(getSnps(mycon))

getSubjects <- function (dbcon) {
    sql <- "SELECT * FROM subjects"
    return(dbGetQuery(dbcon, sql))
}
head(getSubjects(mycon))
class(getSubjects(mycon))
dim(getSubjects(mycon))



########  open the human genome database
humdb <- list.files(system.file("extdata",package="org.Hs.eg.db"),full.names=T)

drv <- SQLite()
con <- dbConnect(drv, dbname=humdb)
dbListTables(con)
dbListFields(con,"genes")
dbListFields(con,"kegg")

dbGetQuery(mycon, sprintf("ATTACH '%s' AS humdb",humdb))

sql <- "SELECT * FROM humdb.genes LIMIT 5"
dbGetQuery(mycon, sql)
sql <- "SELECT * FROM humdb.kegg LIMIT 5"
dbGetQuery(mycon, sql)

sql <- "SELECT * FROM snps,humdb.genes WHERE snps.gene_id == humdb.genes.gene_id LIMIT 5"
temp <- dbGetQuery(mycon, sql)

sql <- "SELECT * FROM snps,humdb.genes,humdb.kegg WHERE snps.gene_id == humdb.genes.gene_id AND humdb.genes._id=humdb.kegg._id"
myanswer <- dbGetQuery(mycon, sql)

#### it looks like it works, but I'm not sure I got the right answer. 

#########  session 4: NetCDF files

######

getGWAScols - three args - nc file, start col, end col

snpfile <- system.file("extdata","snpData.csv",package="AdvancedR2011Data")

system.time(dat601 <- scan(snpfile,what=integer(0),sep=",",skip=600,nlines=1))

snpfilenc <- system.file("extdata","snpData.nc",package="AdvancedR2011Data")
mync <- open.ncdf(snpfilenc)

system.time(dat601nc <- get.var.ncdf( mync , "snpData" , start=c(601,1),count=c(1,-1)))

all.equal(as.numeric(dat601),as.vector(dat601nc),check.attributes=FALSE)

all.equal(as.numeric(dat601),as.vector(dat601nc))

all.equal(dat601,as.vector(dat601nc))

identical(dat601,as.vector(dat601nc))

close.ncdf(mync)


getNcdfVarSummary(mync)


smallsnpfile <- system.file("extdata","small_snpData.nc",package="StudentGWAS")
mync <- open.ncdf(smallsnpfile)

getNcdfVarSummary(mync)



mygetGWAScols <- function (nc,first,last) {
    ### check inputs are of the right class
    if (class(nc) != "ncdf") {
        print (paste(deparse(substitute(nc)),"is not a ncdf object"))
        return(NA)
    }
    if (class(first) != "numeric") {
        print (paste(deparse(substitute(first)),"is not numeric"))
        return(NA)
    }
    if (class(last) != "numeric") {
        print (paste(deparse(substitute(last)),"is not numeric"))
        return(NA)
    }

    #get info on the ncdf for later
    myinfo <- getNcdfVarSummary(nc)

    # check first and last are within bounds and last is bigger than first
    if ((first > myinfo$dim2) | (last > myinfo$dim2) ) {
        print(paste("you're asking for columns",first,"to",last ,"but there are only", myinfo$dim2, "columns"))
        return(NA)
    }
    if (first >= last) {
        print (paste("last must be greater than first. first",first,"last",last))
        return(NA)
    }
    
    dat <- get.var.ncdf( nc , myinfo$name , start=c(1,first),count=c(myinfo$dim1,(last+1-first) ) )

    ###return(as.raw(dat))
    return(dat)
}

dat <- mygetGWAScols(thisnc,1,3)



#########


helpfile <- "/Users/jayoung/temp2/StudentGWAS/man/SQLiteFunctions.Rd"

tools::Rd2HTML (helpfile,out="temp.html")
tools::Rd2txt (helpfile,out="temp.txt")
tools::Rd2latex (helpfile,out="temp.latex")
tools::Rd2ex (helpfile,out="temp.ex")

#########






snps <- matrix(sample((1:3), replace=TRUE, 400), nrow=100, ncol=4) 
nsnp <- ncol(snps) 
nsub <- nrow(snps) 
width <- 2 
delta <- rep.int(0, (nsnp-width)*width)

############## this runs the C routine

out <- .C("composite_linkage_disequilibrium", 
          snp = as.raw(snps),
          n_ind = as.integer(nsub), 
          n_snp = as.integer(nsnp), 
          width = as.integer(width), 
          delta = as.double(delta))
out$delta




.cld <- function(data, width = 5)
{
    nsub <- nrow(data)
    nsnp <- ncol(data)
    if (width > nsnp)
        stop("Width must be less than the number of snps.")
    delta <- rep.int(0, (nsnp-width)*width)
    res <- .C("composite_linkage_disequilibrium",
              snp = as.raw(data),
              n_sub = as.integer(nsub),
              n_snp = as.integer(nsnp),
              width = as.integer(width),
              delta = as.double(delta), PACKAGE="StudentGWAS")
    matrix(res$delta, nrow=(nsnp-width), ncol=width, byrow=TRUE,
           dimnames = list(colnames(data)[1:(nsnp-width)], NULL))
}


x <- as.double(c(1.2,1.44555, 2^0.5))
y <- as.numeric(c(1.2,1.44555, 2^0.5))