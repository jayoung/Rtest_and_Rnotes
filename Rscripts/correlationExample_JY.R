##################
# A = a library pool of 100 clones with uneven coverage (five different coverage values)
A <- rep( c(400,450,500,550,600), 20) 

# b = a library pool of 100 clones with perfectly even coverage
B <- rep( 500, 100  ) 

## pretend we sequence each of those pools twice adding simulated noise (similar noise for each)

# add noise sampled from a normal distribution, with mean 0 and s.d. 25
A_noise1 <- A + rnorm(100, sd=25)
A_noise2 <- A + rnorm(100, sd=25)
B_noise1 <- B + rnorm(100, sd=25)
B_noise2 <- B + rnorm(100, sd=25)

# the replicated samples of library A show very good correlation (0.91 for this sampling)
cor(A_noise1,A_noise2,method="spearman")

# the replicated samples of library B show very poor correlation (0.11 for this sampling), but it doesn't mean this experiment is any worse than the other experiment
cor(B_noise1,B_noise2,method="spearman")

par(mfrow=c(1,2), pty="s") # pty="s" makes plots square
plot( A_noise1,A_noise2, xlim=c(300,700), ylim=c(300,700), pch=19, cex=0.25)
plot( B_noise1,B_noise2, xlim=c(300,700), ylim=c(300,700), pch=19, cex=0.25)
##################


ggpairs(flea, columns = 2:4, 
        lower=list(continuous=wrap("points", 
                                   alpha = 0.1,
                                   xlim=c(100,250))))
        
  
