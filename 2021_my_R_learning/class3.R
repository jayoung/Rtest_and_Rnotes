### lesson 3 tidyverse
library("tidyverse")

## read_csv has an UNDERSCORE!
clinical <- read_csv("data/clinical.csv")

# it's a tibble
clinical 
# displays better than a data.frame
spec(clinical)
# shows class of each column

####### tidyverse functions: column heads don't need quotes

####### select:  takes only some columns
# selecting columns with tidyverse (dplyr)
sel_columns <- select(clinical, tumor_stage, ethnicity, disease)
# in pipe notation:
sel_columns2 <- clinical |> select(tumor_stage, ethnicity, disease)

# range of columns
sel_columns2 <- select(clinical, tumor_stage:vital_status)
sel_columns2a <- select(clinical, 2:4)
sel_columns2b <- select(clinical, starts_with("year"))


##### filter: takes only some rows
filtered_rows <- filter(clinical, disease == "LUSC") 

## filter rows and columns
race_BRCA2 <- select(filter(clinical, disease == "BRCA"), race, ethnicity, disease)

# same task as above, but with pipes
piped <- clinical |>
  select(race, ethnicity, disease) |>
  filter(disease == "BRCA")

temp <- clinical |>
  filter(vital_status=="alive" & cigarettes_per_day<1) |>
  select(gender, years_smoked, year_of_birth)


###### mutate - adding/changing columns
clinical_years <- clinical |>
  mutate(years_to_death = days_to_death / 365)
# transmute is related - drops the existing columns and only keeps the new columns.
clinical_years2 <- clinical_years |>
  transmute(years_to_death = days_to_death / 365)

clinical_years3 <- clinical |>
  mutate(years_to_death = days_to_death / 365,
         months_to_death = days_to_death / 30) |>
  glimpse() # preview data output


temp <- clinical |> 
  filter(disease=="LUSC") |>
  mutate(total_cig=years_smoked*cigarettes_per_day*365) |>
  glimpse()



##### grouping, tallying, summarizing (like tapply):

temp <- clinical |>
  group_by(gender) 

## tally and count are very similar. tally assumes you've already done the grouping
clinical |>
  group_by(gender) |>
  tally() 

clinical |>
  count(gender) 

clinical |>
  count(race, gender) 

# summarize average days to death by gender
clinical |>
  group_by(gender) |>
  summarize(mean_days_to_death = mean(days_to_death, na.rm = TRUE))

smoke_complete <- clinical |>
  filter(!is.na(cigarettes_per_day) & !is.na(age_at_diagnosis))
write_csv(smoke_complete, "data/smoke_complete.csv")

clinical |>
  count(missingCig=is.na(cigarettes_per_day), missingAge=is.na(age_at_diagnosis)) 

clinical |> count(vital_status)

# make sure ALL missing data is removed!
birth_complete <- clinical |>
  filter(!is.na(year_of_birth)) |>
  filter(!is.na(vital_status)) |>
  filter(vital_status != "not reported")

### arrange sorts by column "n"
cancer_counts <- clinical |>
  count(disease) |>
  arrange(n) 

frequent_cancers <- cancer_counts |>
  filter(n >= 500) 

birth_reduced <- birth_complete |>
  filter(disease %in% frequent_cancers$disease)

## write_csv is similar to write.csv but subtly different. No row names. No quotes.
write_csv(birth_reduced, "data/birth_reduced.csv")


clinical |> 
  filter(!is.na(tumor_stage)) |>
  filter(tumor_stage != "not reported") |>
  count(tumor_stage) |>
  arrange(n) |>
  filter(n>200)

clinical |> 
  filter(ethnicity=="hispanic or latino") |>
  filter(race!="white") |>
  count(race)


clinical <- clinical |> mutate(age_at_death=year_of_death-year_of_birth)

clinical |> select(year_of_birth, year_of_death, age_at_death)

## frequent_cancers$disease is a really useful notation - for tibbles it is better than frequent_cancers[,"disease"], which returns another tibble rather than just the data

