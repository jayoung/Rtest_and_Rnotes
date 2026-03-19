######  NGS_functions.R : for functions related to NGS files 

### source this file as follows:
# malik_h_dir <- "/fh/fast/malik_h/"
# if (Sys.info()[["sysname"]]=="Darwin") { malik_h_dir <- "/Volumes/malik_h/" }
# source( paste0(malik_h_dir, "user/jayoung/git_more_repos/Rtest_and_Rnotes/useful_functions/NGS_functions.R") )


########### readFlagstatsFile and combineFlagstatsFiles functions: takes the flagstats file(s) output by samtools flagstats, reads them in, combines multiple samples if relevant, etc.

readFlagstatsFile <- function (myFile, asPercent="no") {
  dat <- scan(myFile, sep="\n", what="character", quiet=TRUE)
  dat <- strsplit(dat, " ")
  datNumQCpass <- as.numeric(sapply(dat,"[[",1))
  if (asPercent=="yes") {
    datNumQCpass <- round(100*datNumQCpass/datNumQCpass[1],1)
  }
  datNumQCfail <- as.numeric(sapply(dat,"[[",3))
  datRest <- lapply( dat, function(x) { x[4:length(x)] } )
  datRest <- sapply( datRest, function(x) { paste(x, collapse=" ") } )
  datRest <- gsub(" \\(QC-passed reads \\+ QC-failed reads\\)","",datRest)
  datRest <- gsub(" \\(.+?nan\\%\\)","",datRest, perl=TRUE)
  datRest <- gsub(" \\(\\d+\\.\\d+\\% : N\\/A\\)","",datRest, perl=TRUE)
  datRest <- gsub(" \\(N\\/A\\ : N\\/A\\)","",datRest, perl=TRUE)
  dat <- data.frame(countQCpass=datNumQCpass, 
                    countQCfail=datNumQCfail, 
                    description=datRest)
  # weirdly, I get the occasional flagstat file that has three extra lines, that all include 'primary' Did I use a different version of samtools, or something? I'm going to drop them
  dat <- dat %>% 
    filter(!str_detect(description, "primary"))
  dat
}

### in Sept 2019 I changed the way I calculate totalSequenced: I now subtract 'supplementary' mappings AS WELL as 'secondary' (previously only subtracted secondary)
# outputFormat can be "original" or "transposed" 
combineFlagstatsFiles <- function (myFiles, 
                                   removeDirFromSampleName=TRUE,
                                   myStripOffTextFromSampleName="", 
                                   asPercent="no",
                                   outputFormat="original") {
  dat <- lapply (myFiles, readFlagstatsFile, asPercent=asPercent)
  names(dat) <- myFiles
  if(removeDirFromSampleName) {
    names(dat) <- sapply( strsplit(names(dat), "/"), function(x) {x[length(x)]} )
  }
  if (myStripOffTextFromSampleName[1] != "") {
    for (thisText in myStripOffTextFromSampleName) {
      names(dat) <- gsub(thisText,"",names(dat))
    }
  }
  checkQCfailReads <- sapply(dat, function(x) { sum(x[,"countQCfail"]) })
  if (sum( checkQCfailReads>0)>0) {
    return("Error! Not expecting to find any reads that failed QC - check everything is OK")
  }
  
  #### go through each sample and add to data.frame
  newDat <- data.frame(description=dat[[1]][,"description"])
  for( x in names(dat)) {
    newDat[,x] <- dat[[x]][,"countQCpass"]
  } 
  rownames(newDat) <- newDat[,"description"]
  
  ### notMapped = in total - mapped  (but I think mapped includes secondaries)
  newDat["notMapped",] <- NA
  newDat["notMapped",2:dim(newDat)[2]] <- newDat["in total",2:dim(newDat)[2]] - newDat["mapped",2:dim(newDat)[2]]
  
  ### totalSequenced = in total - (secondary+supplementary)
  newDat["totalSequenced",] <- NA
  newDat["totalSequenced",2:dim(newDat)[2]] <- newDat["in total",2:dim(newDat)[2]] - (newDat["secondary",2:dim(newDat)[2]] + newDat["supplementary",2:dim(newDat)[2]])
  
  newDat["numReadsMappedAtLeastOnce",] <- NA
  newDat["numReadsMappedAtLeastOnce",2:dim(newDat)[2]] <- newDat["totalSequenced",2:dim(newDat)[2]] - newDat["notMapped",2:dim(newDat)[2]]
  if (length(myFiles) > 1) {
    newDat <- newDat[,2:dim(newDat)[2]]
  } else {
    newColname <- colnames(newDat)[2]
    newDat <- data.frame( myColumn=newDat[,2], row.names=rownames(newDat) )
    colnames(newDat) <- newColname
  }
  if(outputFormat=="transposed") { newDat <- as.data.frame(t(newDat))}
  newDat <- newDat %>% 
    as_tibble(rownames="statistic") %>% 
    pivot_longer(-statistic, names_to="sample", values_to="num_reads")
  return(newDat)
}
