### lesson 2 basic R

dir.create("data")
download.file("https://raw.githubusercontent.com/fredhutchio/R_intro/master/extra/clinical.csv", "data/clinical.csv")

clinical <- read.csv("data/clinical.csv", stringsAsFactors=TRUE)

par(mfrow=c(2,3))
for (i in 1:6) {
  plot(clinical[,i], main=colnames(clinical)[i])
}
