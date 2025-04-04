
#### blank_plot and blank_ggplot - can be useful as placeholders when using patchwork to combine plots

blank_plot <- function() {
    plot(1:10,1:10,"n",bty="n",xaxt="n",yaxt="n", xlab="", ylab="")
}

blank_ggplot <- function() {
    tibble(x=1:10, y=1:10) %>% 
        ggplot() +
        theme_void()
}
