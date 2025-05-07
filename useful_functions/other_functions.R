
#### blank_plot and blank_ggplot - can be useful as placeholders when using patchwork to combine plots

blank_plot <- function() {
    plot(1:10,1:10,"n",bty="n",xaxt="n",yaxt="n", xlab="", ylab="")
}

blank_ggplot <- function() {
    tibble(x=1:10, y=1:10) %>% 
        ggplot() +
        theme_void()
}



### lists the biggest num (default 10) items in memory
# all: can list hidden objects, too
# see also memory.profile()
lsmem <- function (num=10, all=FALSE, units="Kb") {
    z <- sapply(ls(.GlobalEnv), function(x) {
        z1 <- object.size(get(x, envir = .GlobalEnv))
        z1
    })
    z <- as.matrix(rev(sort(z))[1:num])
    didIformat <- "no"
    if (units=="Gb") {z <- round(z/ 1024^3,1) ; didIformat <- "yes"}
    if (units=="Mb") {z <- round(z/ 1024^2,1) ; didIformat <- "yes" }
    if (units=="Kb") {z <- round(z/ 1024,1)   ; didIformat <- "yes" }
    if (units=="bytes") {didIformat <- "yes" }
    if (didIformat=="no") {
        cat("\n\nWARNING - did not recognize the units you asked for. Options are Gb, Mb, Kb\n\n")
        units <- "bytes"
    }
    colnames(z)[1] <- units
    z
}