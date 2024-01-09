### lesson 4 ggplot

library(tidyverse)
birth_reduced <- read_csv("data/birth_reduced.csv")
smoke_complete <- read_csv("data/smoke_complete.csv")

plot(x=smoke_complete$age_at_diagnosis, y=smoke_complete$cigarettes_per_day)

ggplot(data = smoke_complete,
       mapping = aes(x = age_at_diagnosis, y = cigarettes_per_day)) + 
  geom_point() # add a layer of geometry

## same as
ggplot(smoke_complete) + 
  geom_point(aes(x = age_at_diagnosis, y = cigarettes_per_day)) 

ggplot(smoke_complete) +
  geom_point(aes(x = age_at_diagnosis, y = cigarettes_per_day), alpha = 0.1)

ggplot(smoke_complete) +
  geom_point(aes(x = age_at_diagnosis, y = cigarettes_per_day, 
                 color = disease), 
             alpha = 0.1)

ggplot(smoke_complete) +
  geom_point(aes(x = age_at_diagnosis, y = cigarettes_per_day, color = disease), alpha = 0.1) +
  theme_bw() # change background theme

ggplot(smoke_complete) +
  geom_point(aes(x = age_at_diagnosis, y = cigarettes_per_day, color = disease), alpha = 0.1) +
  labs(title = "Age at diagnosis vs cigarettes per day", # title
       x="age (days)", # x axis label
       y="cigarettes per day") +# y axis label
  theme_bw()


ggplot(smoke_complete) +
  geom_point(aes(x = age_at_diagnosis, y = cigarettes_per_day, color = disease)) +
  theme(text = element_text(size = 16)) # increase all font size

#####
ggplot(smoke_complete) +
  geom_point(aes(x = age_at_diagnosis, y = cigarettes_per_day, color = disease), alpha = 0.1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) # rotate and adjust x axis text

ggsave("figures/awesomePlot.jpg", width = 10, height = 10, dpi = 300)

ggsave("figures/awesomePlot.pdf", device="pdf", width = 10, height = 10, dpi = 300)


#####
ggplot(smoke_complete) +
  geom_point(aes(x = age_at_diagnosis, y = years_smoked, color = gender), alpha = 0.1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) # rotate and adjust x axis text

ggsave("figures/awesomePlot.jpg", width = 10, height = 10, dpi = 300)



my_plot <- ggplot(smoke_complete, aes(x = vital_status, y = cigarettes_per_day)) 
my_plot +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.2, color = "purple")


my_plot +
  geom_jitter(alpha = 0.2, color = "purple") + 
  geom_boxplot(outlier.shape = NA) 



yearly_counts <- birth_reduced %>%
  count(year_of_birth, disease) 

ggplot(yearly_counts) +
  geom_line(aes(x = year_of_birth, y = n, 
                color=disease))



yearly_counts_gender <- birth_reduced %>%
  count(year_of_birth, gender) 

ggplot(yearly_counts_gender) +
  geom_line(aes(x = year_of_birth, y = n, 
                color=gender))

### same, but slightly more elegant
birth_reduced %>%
  count(year_of_birth, gender) %>%
  ggplot() +
  geom_line(aes(x = year_of_birth, y = n, 
                lty=gender))


##### faceting.  facet_wrap is low level, facet_grid allows control of rows/columns of plots

ggplot(smoke_complete) +
  geom_point(aes(x = age_at_diagnosis, y = cigarettes_per_day, color = disease)) +
  facet_wrap(vars(disease)) # wraps panels to make a square/rectangular plot

ggplot(smoke_complete) +
  geom_point(aes(x = age_at_diagnosis, y = cigarettes_per_day, color = disease)) +
  facet_wrap(vars(disease)) # wraps panels to make a square/rectangular plot


# facet_grid
ggplot(smoke_complete) +
  geom_point(aes(x = age_at_diagnosis, y = cigarettes_per_day, color = disease)) +
  facet_grid(rows = vars(vital_status), cols=vars(disease)) 

birth_reduced %>%
  count(year_of_birth, gender) %>%
  ggplot() +
  geom_line(aes(x = year_of_birth, y = n)) +
  facet_wrap(vars(gender), ncol=1)


# ~gender works, as does vars(gender) and "gender".  plain gender does not. 
birth_reduced %>%
  count(year_of_birth, gender) %>%
  ggplot() +
  geom_line(aes(x = year_of_birth, y = n)) +
  facet_wrap(~gender, ncol=1)

birth_reduced %>%
  count(year_of_birth, gender) %>%
  ggplot() +
  geom_line(aes(x = year_of_birth, y = n)) +
  facet_wrap("gender", ncol=1)
