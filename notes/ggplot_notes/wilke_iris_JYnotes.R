# we got three regression lines because we used color

ggplot(data=iris, aes(Sepal.Length, Sepal.Width, color=Species)) + 
  geom_point() +
  geom_smooth()



# regress on all data
ggplot(data=iris, aes(Sepal.Length, Sepal.Width, color=Species)) + 
  geom_point() +
  geom_smooth(aes(group=1))
