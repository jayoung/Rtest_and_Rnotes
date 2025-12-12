# file manipulation:

# example data:
dir.create("data")
download.file("https://raw.githubusercontent.com/fredhutchio/R_intro/master/extra/clinical.csv", "data/clinical.csv")

#### read.csv (regular R) versus read_csv (tidyverse)
clinical_old <- read.csv("data/clinical.csv", stringsAsFactors=TRUE)
clinical <- read_csv("data/clinical.csv")

# also write.csv versus write_csv

#### viewing objects
str(clinical)
# str : compactly Display the Structure of an Arbitrary R Object

# summary can also be run on a data frame, and it runs on each column

## Rstudio projects:
# https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects

### recommended project organization:
# data_raw folder for input files
# data for output files

## other Rstudio tips:
#  https://garthtarr.github.io/avfs/tips.html
# for example, to get "<-" (which I always make typos in) you can do control and the - key 
# and control-1 and control-2 switch you between the script and the console sections of the window 


###### key tidyverse functions:


# |>    # pipe.  shortcut =command-shift-M
spec()   # on a tibble, shows the class of each column
select() # chooses columns
filter() # chooses rows
mutate() # add extra columns
transmute() # creates new columns and drops the existing ones
glimpse() # can pipe into this to get a look at the end result

rename() # to rename columns
rename_with() # allows column renaming using functions like gsub, toupper, etc

n_distinct() # similar to length(unique())

###### group_by / tally / count / summarize

## tally needs to be run after group_by
clinical |>
  group_by(gender) |>
  tally() 

## count does not need to be run after group_by
clinical |>
  count(gender)

clinical |>
  count(gender, disease)

clinical |>
  count(disease, gender)

## group_by doesn't change the object, but it does changes how objects interact with downstream functions: example
by_cyl <- mtcars |> group_by(cyl)

# grouping doesn't change how the data looks (apart from listing
# how it's grouped):
by_cyl

# It changes how it acts with the other dplyr verbs:
by_cyl |> summarise(
  disp = mean(disp),
  hp = mean(hp)
)
by_cyl |> filter(disp == max(disp))

# compare that with 
mtcars |> summarise(
  disp = mean(disp),
  hp = mean(hp)
)
mtcars |> filter(disp == max(disp))

# Each call to summarise() removes a layer of grouping
by_vs_am <- mtcars |> group_by(vs, am)
by_vs <- by_vs_am |> summarise(n = n())
by_vs
by_vs |> summarise(n = sum(n))

# To removing grouping, use ungroup
by_vs |>
  ungroup() |>
  summarise(n = sum(n))



## summarize (similar to tapply)
clinical |>
  group_by(gender) |>
  summarize(mean_days_to_death = mean(days_to_death, na.rm = TRUE))

clinical |>
  group_by(gender) |>
  summarize(mean_days_to_death = mean(days_to_death, na.rm = TRUE),
            numPeople=n())


####### arrange to sort by a column:
clinical |>
  count(disease)

clinical |>
  count(disease) |>
  arrange(n) 
# arrange(desc(n)) for reverse order

###### in the tidyverse, the $ notation for columns could be really useful for tibbles: it is better than frequent_cancers$disease yields the data as a vector, wherease frequent_cancers[,"disease"], returns a 1-column tibble



###### ggplot

# remember to use alpha=0.1 for transparency in geom_point()

## save current plot as jpg/pdf etc
ggsave("figures/awesomePlot.pdf", device="pdf", width = 10, height = 10, dpi = 300)

## save a plot in memory:
ggsave("figures/awesomePlot.pdf", plot=myplot, device="pdf", width = 10, height = 10, dpi = 300)

## popular themes:
theme_bw()
theme_minimal() 
theme_light() 
theme_classic()
theme_void() # can be useful as a starting point to create a new hand-crafted theme.


### patchwork package - yet another way to put multiple ggplots on the same page
# After you have loaded the patchwork package you can use:
# + to combine plots in a way that patchwork will determine (e.g. 4 plots go in a 2x2 grid)
# | to place plots next to each other, 
# / to arrange them vertically, and plot_layout() to determine how much space each plot uses
# use parentheses to group plots and nest layouts
  
library(patchwork)

plot1 <- ggplot(data = surveys_complete, aes(x = species_id, y = weight)) +
  geom_boxplot() +
  labs(x = "Species", y = expression(log[10](Weight))) +
  scale_y_log10()

plot2 <- ggplot(data = yearly_counts, aes(x = year, y = n, color = genus)) +
  geom_line() + 
  labs(x = "Year", y = "Abundance")

plot1 / plot2 + plot_layout(heights = c(3, 2))

p1 <- ggplot(mtcars) + geom_point(aes(mpg, disp))
p2 <- ggplot(mtcars) + geom_boxplot(aes(gear, disp, group = gear))

p1 + p2

p3 <- ggplot(mtcars) + geom_smooth(aes(disp, qsec))
p4 <- ggplot(mtcars) + geom_bar(aes(carb))

# these two things give the same result:
# 1
(p1 | p2 | p3) /
  p4
# 2
(p1 + p2 + p3) /
  p4

plot_annotation()  # to put things in the 'outer' margins or to label each plot

####### what next?  I want to understand long versus wide datasets in the tidyverse
## https://datacarpentry.org/R-ecology-lesson/03-dplyr.html
# Describe the concept of a wide and a long table format and for which purpose those formats are useful.
# Describe what key-value pairs are.
# Reshape a data frame from long to wide format and back with the spread and gather commands from the tidyr package.




