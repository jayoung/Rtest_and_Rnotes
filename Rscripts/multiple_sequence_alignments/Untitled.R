# Libraries
library(ggplot2)
library(dplyr)
library(patchwork) # To display 2 charts together


# Build dummy data
data <- data.frame(
    day = as.Date("2019-01-01") + 0:99,
    temperature = runif(100) + seq(1,100)^2.5 / 10000,
    price = runif(100) + seq(100,1)^1.5 / 10
)

# Most basic line chart
p1 <- ggplot(data, aes(x=day, y=temperature)) +
    geom_line(color="#69b3a2", linewidth=2) +
    ggtitle("Temperature: range 1-10") 

p2 <- ggplot(data, aes(x=day, y=price)) +
    geom_line(color="grey",linewidth=2) +
    ggtitle("Price: range 1-100") 

# Display both charts side by side thanks to the patchwork package
p1 + p2



### same plot diff axes
# Value used to transform the data
coeff <- 10

ggplot(data, aes(x=day)) +
    
    geom_line( aes(y=temperature)) + 
    geom_line( aes(y=price / coeff)) + # Divide by 10 to get the same range than the temperature
    
    scale_y_continuous(
        
        # Features of the first axis
        name = "First Axis",
        
        # Add a second axis and specify its features
        sec.axis = sec_axis(~.*coeff, name="Second Axis")
    )