library(tidyverse)

####
##### https://datacarpentry.org/R-ecology-lesson/03-dplyr.html
# Describe the concept of a wide and a long table format and for which purpose those formats are useful.
# Describe what key-value pairs are.
# Reshape a data frame from long to wide format and back with the spread and gather commands from the tidyr package.


# dplyr package - for working with dataframes better for large data (can use relational databases rather than storing everythin in memory)
# tidyr - manipulate dataframes (e.g. long versus wide)


### example data: animal species diversity and weights found within plots at our study site.
#download.file(url = "https://ndownloader.figshare.com/files/2292169",
#              destfile = "data/portal_data_joined.csv")

surveys <- read_csv("data/portal_data_joined.csv")
# each row is an observation

### character columns do not produce a useful output with 'summary', but factors do:
summary(surveys$sex)

surveys$sex <- factor(surveys$sex)
summary(surveys$sex)

# Change the columns taxa and genus in the surveys data frame into a factor.
surveys$taxa <- factor(surveys$taxa)
surveys$genus <- factor(surveys$genus)

summary(surveys)

# How many rabbits were observed?  75
surveys |> filter(taxa=="Rabbit") |> nrow

# How many different genera are in the genus column?  26
nlevels(surveys$genus)
surveys |> count(genus) |> nrow


#### converting factors to numeric - be careful
year_fct <- factor(c(1990, 1983, 1977, 1998, 1990))
as.numeric(year_fct)               # Wrong! And there is no warning...
as.numeric(as.character(year_fct)) # Works...
as.numeric(levels(year_fct))[year_fct]    # The recommended way.

#### convenient plotting for factors:
plot(surveys$sex)

# that ignores the 1748 NA individuals 
surveys |> count(sex)
levels(surveys$sex)

# add NA as a level to be considered, and rename it. This is basic R, not tidyr:
sex <- surveys$sex
sex <- addNA(sex)
levels(sex)
# renaming a single level of a factor changes the underlying data, too
levels(sex)[3] <- "undetermined"
levels(sex)

plot(sex)

# Rename “F” and “M” to “female” and “male” respectively.
levels(sex)[1:2] <- c("female","male")

### recreate the barplot such that “undetermined” is first (before “female”)?

head(sex,20)
# [1] male         male         undetermined undetermined undetermined undetermined undetermined
# [8] undetermined male         undetermined undetermined male         male         male        
# [15] male         undetermined male         male         male         male        
# Levels: female male undetermined

table(sex)
# sex
# female         male undetermined 
# 15690        17348         1748 

## suprisingly I do not need to include 'as.character' here - it doesn't screw things up. However, if I were to have missed out a category in my list of levels it would drop those data.points
sex2 <- factor(sex, levels=c("undetermined","female", "male"))
table(sex2)
#undetermined       female         male 
#        1748        15690        17348 
plot(sex2)

## relevel - a quick way to put a single level first (often because one level is the 'reference' level, but not in this case):
sex3 <- relevel(sex, "undetermined")
table(sex3)

## To quickly reverse the order of all levels in a factor:
sex4 <- factor(sex, levels=rev(levels(sex)))
table(sex4)


###### date info:
# the best practice for dealing with date data is to ensure that each component of your date is stored as a separate variable
str(surveys)
# we have month, day, year columns - that's good

# a nice package for date manipulation
library(lubridate)

# lubridate's ymd() function takes a vector representing year, month, and day, and converts it to a Date vector. Date is a class of data recognized by R as being a date and can be manipulated as such. The argument that the function requires is flexible, but, as a best practice, is a character vector formatted as “YYYY-MM-DD”.

# Let’s create a date object and inspect the structure:
  
my_date <- ymd("2015-01-01")
str(my_date)

surveys$date <- ymd(paste(surveys$year, surveys$month, surveys$day, sep = "-"))
summary(surveys$date)
#         Min.      1st Qu.       Median         Mean      3rd Qu.         Max.         NA's 
# "1977-07-16" "1984-03-12" "1990-07-22" "1990-12-15" "1997-07-29" "2002-12-31"        "129" 

# 129 rows failed to parse - why?  fix it.
missing_dates <- surveys[is.na(surveys$date), c("year", "month", "day")]
head(missing_dates)

##  all are day 31, some are month 4, some month 9. Both those months have 30 days.  Fix the days!
table(missing_dates[["day"]])
#  31 
# 129 
table(missing_dates[["month"]])
#  4  9 
# 70 59 

missing_dates <- missing_dates |> mutate(day_fixed=day-1)
missing_dates <- missing_dates |> 
  mutate(date=ymd(paste(year, month, day_fixed, sep = "-")  ))
# add fixed dates back to original object:
surveys[is.na(surveys$date), "date"] <- missing_dates$date


##### remember the negative notation to remove some columns
select(surveys, -record_id, -species_id) 


##### Challenge. Create a new data frame from the surveys data that meets the following criteria: 
# contains only the species_id column 
# and a new column called hindfoot_cm containing the hindfoot_length values converted to centimeters (it's in mm now)
# In this hindfoot_cm column, there are no NAs and all values are less than 3.
# Hint: think about how the commands should be ordered to produce this data frame!

surveys2 <- surveys |> 
  mutate(hindfoot_cm=hindfoot_length/10) |> 
  filter(hindfoot_cm<3) |> 
  filter(!is.na(hindfoot_cm)) |> 
  select(species_id)

###### grouping, summarizing
surveys |>
  group_by(sex) |>
  summarize(mean_weight = mean(weight, na.rm = TRUE))

# there are NAs at the bottom
surveys |>
  group_by(sex, species_id) |>
  summarize(mean_weight = mean(weight, na.rm = TRUE)) |> 
  tail()

# get >1 summary statistic for each group add arrange to sort
surveys |>
  filter(!is.na(weight)) |>
  filter(!is.na(sex)) |>
  group_by(sex, species_id) |>
  summarize(mean_weight = mean(weight),
            min_weight = min(weight)) |> 
  arrange(min_weight) 

# descending sort
surveys |>
  filter(!is.na(weight)) |>
  filter(!is.na(sex)) |>
  group_by(sex, species_id) |>
  summarize(mean_weight = mean(weight),
            min_weight = min(weight)) |> 
  arrange(desc(min_weight) )


# these are the same:
surveys |>
  count(sex) 
  
surveys |>
  group_by(sex) |>
  summarise(count = n())

# alternative to arrange:
surveys |>
  count(sex, sort=TRUE) 

# sorting on two columns:
surveys |>
  count(sex, species) |>
  arrange(species, desc(n))


### challenge:  Use group_by() and summarize() to find the mean, min, and max hindfoot length for each species (using species_id). Also add the number of observations (hint: see ?n).
surveys |> 
  filter(!is.na(hindfoot_length)) |> 
  group_by(species_id) |> 
  summarise(meanFoot=mean(hindfoot_length),
            minFoot=min(hindfoot_length),
            maxFoot=max(hindfoot_length),
            numberAnimals=n())


### challenge:  What was the heaviest animal measured in each year? Return the columns year, genus, species_id, and weight.  This is almost it, although if there is a tie in a year I keep both heaviest animals
surveys |> 
  group_by(year) |> 
  filter(!is.na(weight)) |> 
  filter(weight==max(weight)) |> 
  select(year, genus, species_id, weight) |> 
  arrange(year)


surveys_gw <- surveys |> 
  filter(!is.na(weight)) |> 
  group_by(plot_id, genus) |> 
  summarise(meanWt=mean(weight, na.rm=TRUE))

### spread turns a LONG tibble into a WIDE tibble
# somehow it knew to keep the first column (plot_id) the same.
# then the 'key' argument tells it to take the second column (genus) and make those new column headers, using meanWt as the values in those columns
surveys_spread <- surveys_gw |>
  spread(key = genus, value = meanWt)

### turns out spread has been superceded by pivot_wider: this gives the SAME result
surveys_spread2 <- surveys_gw |>
  pivot_wider(names_from = genus, values_from = meanWt)
# using values_fill we can turn the NAs into 0
surveys_spread2 <- surveys_gw |>
  pivot_wider(names_from = genus, values_from = meanWt, values_fill=0)

### gather does the opposite -turns a WIDE tibble into a LONG tibble 
surveys_gather <- surveys_spread |>
  gather(key = "genus", value = "mean_weight", -plot_id)
# -plot_id means take every column except plot_id and treat it as a genus value and put the values in a mean_weight column.  It does keep the plot_id as the first column, just like spread did.

### and gather has been superceded by pivot_longer. This gives identical results
surveys_gather2 <- surveys_spread |> 
  pivot_longer(cols=-plot_id, names_to = "genus", values_to = "mean_weight") |> 
  arrange(genus)
# maybe we want to drop the NAs to resemble the original summarized data, or maybe we don't.  a quick spread-gather can be useful if you want to ensure all combinations are present and create the relevant NA values
surveys_gather2 <- surveys_spread |> 
  pivot_longer(cols=-plot_id, names_to = "genus", values_to = "mean_weight", values_drop_na=TRUE) |> 
  arrange(genus)


#### challenges:
### 1. Spread the surveys data frame with:
# year as columns, plot_id as rows, and the number of genera per plot as the values. 
# You will need to summarize before reshaping, and use the function n_distinct() to get the number of unique genera within a particular chunk of data. It’s a powerful function! See ?n_distinct for more

# num genera per plot per year
surveys_spread_genera <- surveys |> 
  group_by(year, plot_id) |> 
  summarise(numGenera=n_distinct(genus)) |> 
  pivot_wider(id_cols=plot_id, names_from=year, values_from=numGenera) 

# works - I did some spot checks on that, e.g.
temp <- surveys |> 
  filter(year==1984, plot_id==3) 
n_distinct( temp$genus)
rm(temp)

### 2. Now take that data frame and gather() it again, so each row is a unique plot_id by year combination.
surveys_spread_genera |> 
  pivot_longer(cols=-plot_id, names_to="year", values_to="numGenera") |> 
  arrange(year, plot_id) |> 
  select(year, plot_id, numGenera) # select is just to reorder columns so it looks more like the output from the summarise  part of step 1 of this challenge

### 3. The surveys data set has two measurement columns: hindfoot_length and weight. 
# This makes it difficult to do things like look at the relationship between mean values of each measurement per year in different plot types. 
# Let’s walk through a common solution for this type of problem. 
# First, use gather() to create a dataset where we have a key column called measurement and a value column that takes on the value of either hindfoot_length or weight. 
# Hint: You’ll need to specify which columns are being gathered.

surveys_tricky <- surveys |> 
  #select(year,plot_type,hindfoot_length,weight)
  pivot_longer(cols=c(hindfoot_length,weight),names_to="measurement") |> 
  arrange(measurement) |> # cosmetic, so I get the identical answer to the website
  glimpse()
  
# the answer from the website:
surveys_long <- surveys |>
  gather("measurement", "value", hindfoot_length, weight)

identical(surveys_tricky, surveys_long)
# TRUE !

### 4. With this new data set, calculate the average of each measurement in each year for each different plot_type. 
# Then spread() them into a data set with a column for hindfoot_length and weight. 
# Hint: You only need to specify the key and value columns for spread().

temp1 <- surveys_tricky |> 
  group_by(year, plot_type, measurement) |> 
  summarise(mean=mean(value, na.rm=TRUE)) |> 
  pivot_wider(id_cols=-mean, names_from=measurement, values_from=mean)

# answer from the website:
temp2 <- surveys_long |>
  group_by(year, measurement, plot_type) |>
  summarize(mean_value = mean(value, na.rm=TRUE)) |>
  spread(measurement, mean_value)

identical(temp1, temp2)
# FALSE   ## the Groups component looks different:
# temp1:  # Groups:   year, plot_type [130]
# temp2:  # Groups:   year [26]
identical(as.data.frame(temp1), as.data.frame(temp2))
# TRUE


### clean up missing data
surveys_complete <- surveys |>
  filter(!is.na(weight),           # remove missing weight
         !is.na(hindfoot_length),  # remove missing hindfoot_length
         !is.na(sex))                # remove missing sex

#  remove observations for rare species (i.e., that have been observed less than 50 times). We will do this in two steps: first we are going to create a data set that counts how often each species has been observed, and filter out the rare species; then, we will extract only the observations for these more common species:

## Extract the most common species_id
species_counts <- surveys_complete |>
  count(species_id) |> 
  filter(n >= 50)

## Only keep the most common species
surveys_complete <- surveys_complete |>
  filter(species_id %in% species_counts$species_id)

# only rodents survived the filtering:
surveys_complete |> count(taxa)
surveys_complete |> count(species)
surveys_complete |> count(genus)

########### plotting

surveys_plot <- ggplot(data=surveys_complete, 
       aes(x = weight, y = hindfoot_length)) 

surveys_plot + geom_point()
surveys_plot + geom_point(aes(col=genus))

surveys_plot + geom_point(alpha=0.1)

library(hexbin)
surveys_plot + geom_hex()

## Use what you just learned to create a scatter plot of weight over species_id with the plot types showing in different colors. Is this a good way to show this type of data?
ggplot(data=surveys_complete, 
       aes(x = species_id , y = weight)) +
  geom_point(aes(col=plot_type))

ggplot(data = surveys_complete, mapping = aes(x = species_id, y = weight)) +
  geom_boxplot()

ggplot(data = surveys_complete, mapping = aes(x = species_id, y = weight)) +
  geom_jitter(alpha = 0.3, color = "tomato") +
  geom_boxplot(alpha = 0) 


### Challenges
# Boxplots are useful summaries, but hide the shape of the distribution. For example, if there is a bimodal distribution, it would not be observed with a boxplot. An alternative to the boxplot is the violin plot (sometimes known as a beanplot), where the shape (of the density of points) is drawn.
# Replace the box plot with a violin plot; see geom_violin().

ggplot(data = surveys_complete, mapping = aes(x = species_id, y = weight)) +
  geom_violin(scale="width")


# In many types of data, it is important to consider the scale of the observations. For example, it may be worth changing the scale of the axis to better distribute the observations in the space of the plot. Changing the scale of the axes is done similarly to adding/modifying other components (i.e., by incrementally adding commands). Try making these modifications:
#   Represent weight on the log10 scale; see scale_y_log10().

ggplot(data = surveys_complete, mapping = aes(x = species_id, y = weight)) +
  geom_violin(scale="width") +
  scale_y_log10()

# So far, we’ve looked at the distribution of weight within species. Try making a new plot to explore the distribution of another variable within each species.
ggplot(data = surveys_complete, mapping = aes(x = species_id)) +
  geom_bar(stat="count", aes(fill=sex))

# Create boxplot for hindfoot_length. Overlay the boxplot layer on a jitter layer to show actual measurements.
# Add color to the data points on your boxplot according to the plot from which the sample was taken (plot_id). 
ggplot(data = surveys_complete, mapping = aes(x = species_id, y=hindfoot_length)) +
  geom_jitter(alpha=0.1, aes(color=factor(plot_id))) +
  geom_boxplot()

yearly_counts <- surveys_complete |>
  count(year, genus)

ggplot(data = yearly_counts, aes(x = year, y = n, color=genus)) +
  geom_line()

# or, using pipe:
yearly_counts |> 
  ggplot(mapping = aes(x = year, y = n, color = genus)) +
  geom_line()

## facets:
ggplot(data = yearly_counts, aes(x = year, y = n)) +
  geom_line() +
  facet_wrap(facets = vars(genus))

# Now we would like to split the line in each plot by the sex of each individual measured. 
surveys_complete |> 
  count(year, genus, sex) |> 
  ggplot(aes(x = year, y = n, colour=sex)) +
    geom_line() +
    facet_wrap(facets = vars(genus)) +
    theme_void()

## Challenge: Use what you just learned to create a plot that depicts how the average weight of each species changes through the years.
surveys_complete |> 
  group_by(species,year) |> 
  summarize(meanWt=mean(weight)) |> 
  ggplot(aes(x=year, y=meanWt, color=species)) +
  geom_line() +
  facet_wrap(vars(species))

### cosmetics:
yearly_sex_counts <- surveys_complete |>
  count(year, genus, sex)
ggplot(data = yearly_sex_counts, mapping = aes(x = year, y = n, color = sex)) +
  geom_line() +
  facet_wrap(vars(genus)) +
  labs(title = "Observed genera through time",
       x = "Year of observation",
       y = "Number of individuals") +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, 
                                   hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        strip.text = element_text(face = "italic"),
        text = element_text(size = 16))

### or, similar thing using a custom theme:
grey_theme <- theme(axis.text.x = element_text(colour="grey20", size = 12, 
                                               angle = 90, hjust = 0.5, 
                                               vjust = 0.5),
                    axis.text.y = element_text(colour = "grey20", size = 12),
                    text=element_text(size = 16))

ggplot(surveys_complete, aes(x = species_id, y = hindfoot_length)) +
  geom_boxplot() +
  grey_theme

ggplot(data = yearly_sex_counts, mapping = aes(x = year, y = n, color = sex)) +
  geom_line() +
  facet_wrap(vars(genus)) +
  labs(title = "Observed genera through time",
       x = "Year of observation",
       y = "Number of individuals") +
  theme_bw() + 
  grey_theme




##### merging/joining datasets
# https://m-clark.github.io/data-processing-and-visualization/tidyverse.html#merging-data
# inner_join: return all rows from x where there are matching values in y, and all columns from x and y. If there are multiple matches between x and y, all combination of the matches are returned.
# left_join: return all rows from x, and all columns from x and y. Rows in x with no match in y will have NA values in the new columns. If there are multiple matches between x and y, all combinations of the matches are returned.
# right_join: return all rows from y, and all columns from x and y. Rows in y with no match in x will have NA values in the new columns. If there are multiple matches between x and y, all combinations of the matches are returned.
# semi_join: return all rows from x where there are matching values in y, keeping just columns from x. It differs from an inner join because an inner join will return one row of x for each matching row of y, where a semi join will never duplicate rows of x.
# anti_join: return all rows from x where there are not matching values in y, keeping just columns from x.
# full_join: return all rows and all columns from both x and y. Where there are not matching values, returns NA for the one missing.

# Probably the most common is a left join, where we have one primary data set, and are adding data from another source to it while retaining it as a base. The following is a simple demonstration.

band_members
band_instruments

left_join(band_members, band_instruments)
right_join(band_members, band_instruments)
inner_join(band_members, band_instruments)
full_join(band_members, band_instruments)
anti_join(band_members, band_instruments)
anti_join(band_instruments, band_members)


