# https://m-clark.github.io/data-processing-and-visualization/tidyverse.html#tidyverse-exercises

# Install and load the dplyr ggplot2movies packages. Look at the help file for the movies data set, which contains data from IMDB.

# install.packages('ggplot2movies')
library(ggplot2movies)
data('movies')

# Exercise 1a
# Use mutate to create a centered version of the rating variable. A centered variable is one whose mean has been subtracted from it. 
movies %>% 
  mutate(rating_norm=rating-mean(rating)) %>% 
  str

# Exercise 1b
# Use filter to create a new data frame that has only movies from the years 2000 and beyond. 
movies %>% 
  mutate(rating_norm=rating-mean(rating)) %>% 
  filter(year>=2000) %>% 
  str

# Exercise 1c
# Use select to create a new data frame that only has the title, year, budget, length, rating and votes variables. There are at least 3 ways to do this.
movies %>% 
  select(title, year, budget, length, rating, votes)

# Exercise 1d
# Rename the length column to length_in_min (i.e. length in minutes).
movies %>% 
  rename(length_in_min=length) 


# Exercise 2
# Use group_by to group the data by year, and summarize to create a new variable that is the average budget. The summarize function works just like mutate in this case.  Instead of doing na.rm=TRUE on the mean column, I get rid of those movies entirely, so the year doesn't show up if there aren't any movies whose budget is known
movies %>% 
  filter(!is.na(budget)) %>% 
  group_by(year) %>% 
  summarise(meanBudget=mean(budget))



# Exercise 3
# Use pivot_longer to create a ‘tidy’ data set from the following.
dat <- tibble(id = 1:10,
             x = rnorm(10),
             y = rnorm(10))
dat_tidy <- dat %>% 
  pivot_longer(cols=-id, names_to="measurement")


# Exercise 4
# Now put several actions together in one set of piped operations.
# Filter movies released after 1990
# select the same variables as before but also the mpaa, Action, and Drama variables
# group by mpaa and (your choice) Action or Drama
# get the average rating
movies %>% 
  filter(year>1990) %>% 
  select(title, year, budget, length, rating, votes, mpaa, Action, Drama) %>% 
  group_by(mpaa,Action) %>% 
  summarise(aveRating=mean(rating)) 


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


###### unite/separate


stocks <- data.frame(
  time = as.Date('2009-01-01') + 0:9,
  X = rnorm(10, 0, 1),
  Y = rnorm(10, 0, 2),
  Z = rnorm(10, 0, 4)
)

stocks %>% head

stocks %>% 
  pivot_longer(
    cols      = -time,   # works similar to using select()
    names_to  = 'stock', # the name of the column that will have column names as labels
    values_to = 'price'  # the name of the column for the values
  ) %>% 
  head()

stocks <- data.frame(
  time = as.Date('2009-01-01') + 0:9,
  X_1 = rnorm(10, 0, 1),
  X_2 = rnorm(10, 0, 1),
  Y_1 = rnorm(10, 0, 2),
  Y_2 = rnorm(10, 0, 2),
  Z_1 = rnorm(10, 0, 4),
  Z_2 = rnorm(10, 0, 4)
)

head(stocks)

stocks %>% 
  pivot_longer(
    cols = -time,
    names_to = c('stock', 'entry'),   ### two names because we're splitting column headers
    names_sep = '_',    ### this does some splitting on the column headers
    values_to = 'price'
  ) %>% 
  head()


#### get bball example data:
library(rvest)
current_year = lubridate::year(Sys.Date())
url = glue::glue("http://www.basketball-reference.com/leagues/NBA_{current_year-1}_totals.html")
bball = read_html(url) %>% 
  html_nodes("#totals_stats") %>% 
  html_table() %>% 
  data.frame() 


#### mutate can be used to change class of columns, not just to add columns:
# across means work on multiple columns
bball = bball %>% 
  mutate(across(c(-Player, -Pos, -Tm), as.numeric))
# I get warnings about NA coercion, because some values were ""

## some more examples
bball = bball %>% 
  mutate(
    trueShooting = PTS / (2 * (FGA + (.44 * FTA))),
    effectiveFG  = (FG + (.5 * X3P)) / FGA,
    shootingDif  = trueShooting - FG.
  )

summary(select(bball, shootingDif))

bball %>% 
  separate(Player, into=c('first_name', 'last_name'), sep=' ') %>% 
  select(1:5) %>% 
  head()
# there were warnings about players who had less or more than two names

## unite combines column values (default sep="_") (it also gets rid of the original columns unless remove=FALSE)
