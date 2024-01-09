library("RSQLite")
library("ncdf")

### definition
setClass("GWASdata", 
     representation(  datapath="character", 
                      dataconn="list",
                      metadatapath="character", 
                      metadataconn="SQLiteConnection",
                      nrow = "integer",
                      ncol = "integer"
                    ) )

showClass("GWASdata")

#### constructor

GWASdata <- function( datapath, metadatapath ) {
     ### open netCDF  
     dataconn <- open.ncdf(datapath)
     class(dataconn) <- NULL
     
     ### open SQLite
     drv <- SQLite()
     metadataconn <- dbConnect(drv, dbname=metadatapath)

     ### count number of subjects and assign to nrow
     sql <- "SELECT id FROM subjects"   ### SELECT count(*) FROM snps - count is a standard SQL function. then just need to do [[1]] on the result, no need to take nrow.  better because doesn't create a big table taking up memory
     nrow <- nrow(dbGetQuery(metadataconn, sql))

     ### count snps in the sql tables and assign to ncol
     sql <- "SELECT id FROM snps"
     ncol <- nrow(dbGetQuery(metadataconn, sql))

     ### now call new, now that we have all those things
     new("GWASdata",  datapath=datapath, 
                      dataconn=dataconn,
                      metadatapath=metadatapath, 
                      metadataconn=metadataconn,
                      nrow=nrow,
                      ncol=ncol)
                      
    ##### should we do error catching here?  maybe so.  see if the funcitons inside give useful error messages - if not I should add my own checks/messages. don't duplicate things that go in the validity method, though
}

### test the constructor
datapath <- system.file("extdata","small_snpData.nc",package="StudentGWAS")

metadatapath <- system.file("extdata","small_metadata.sqlite",package="StudentGWAS")

gwas <- GWASdata(datapath,metadatapath)

##### define accessors

setGeneric("dataconn",function (x) standardGeneric("dataconn") )
setMethod("dataconn","GWASdata", function (x) {
    dataconn <- x@dataconn
    class(dataconn) <- "ncdf"
    return(dataconn)
    } )

setMethod("nrow","GWASdata", function (x) x@nrow)

setMethod("ncol","GWASdata", function (x) x@ncol)

setMethod("dim","GWASdata", function (x) {
    return ( c(nrow(x),ncol(x) ) )
    ### Herve's code didn't have the return. can I do without?  is it because I have the {}?
} )


tempfun <- function (x) { y<-x*x ; y }

tempfun2 <- function (x) y<-x*x ; y


###### test accessors
dataconn(gwas)
nrow(gwas)
ncol(gwas)
dim(gwas)

#######   exercise 4


setGeneric("getCols", signature="x", function (x, first, last=first) standardGeneric("getCols") )

#### signature is the object I'm dispatching on.  Specify it here, because I don't want first and last to be dispatched on as they're options, not objects

setMethod("getCols", "GWASdata", function (x, first, last=first) {
    return ( getGWAScols ( dataconn(x),first,last) )
} )

getCols(gwas,3,4)

############ show method

setMethod("show","GWASdata",function (object) 
    cat(class(object),"instance with",nrow(object),"subjects and",ncol(object),"SNPs")
 )

show(gwas)

#### herve's vaildity method is split using a couple of helper functions - he likes not to let functions get too long. He does a very thorough check on this stuff.  look at his code.

setAs("GWASdata","matrix", function(from) {
        ###do something to the from
        getCols( from,1,ncols(from) )
    }
)