#### chapter 1

library(lattice)
data(Chem97,package="mlmRev")
data(gvhd10,package="latticeExtra")

####### prepare the data for lattice:

## reshape function:

summary(Indometh)
wide <- reshape(Indometh, v.names = "conc", idvar = "Subject",timevar = "time", direction = "wide")
wide

long <- reshape(wide, direction = "long")
long2 <- reshape(wide, idvar = "Subject", varying = list(2:12), v.names = "conc", direction = "long")


################### making tables

xtabs(~score,data=Chem97)
table(Chem97$score)

#two-way table
xtabs(~ score + gender,data=Chem97)
xtabs(~ gender + score,data=Chem97)
table(Chem97$score,Chem97$gender)

#three-way table: ftable makes things look nice.  
#ftable on xtabs looks a little bit nicer than on table, because it has the variable labels
 
xtabs(~ gender + score + age,data=Chem97)
ftable(xtabs(~ gender + score + age,data=Chem97))
ftable(Chem97$gender,Chem97$score,Chem97$age)
ftable(table(Chem97$gender,Chem97$score,Chem97$age))

############ now some lattice plotting examples

### these use score as a conditioning variable (produces separate plots)
histogram(~ gcsescore | factor(score),data=Chem97)

densityplot(~ gcsescore | factor(score),data=Chem97,plot.points=FALSE,ref=TRUE)
### ref adds the 0 line

histogram(~ gcsescore | score,data=Chem97)

### play with plot labels
histogram(~ gcsescore | factor(score),data=Chem97,strip=strip.custom(strip.names=TRUE,var.name=Score) )

### now use score as a grouping variable to plot on same graph
densityplot(~ gcsescore,data=Chem97, groups=score, plot.points=FALSE, ref=TRUE, auto.key=list(columns=3))

histogram(~ gcsescore,data=Chem97, groups=score, plot.points=FALSE, ref=TRUE, auto.key=list(columns=3))

### xyplot - tile different cuts of the data
Depth <- equal.count(quakes$depth, number=8, overlap=.1)
xyplot(lat ~ long | Depth, data = quakes)
update(trellis.last.object(),strip = strip.custom(strip.names = TRUE, strip.levels = TRUE),par.strip.text = list(cex = 0.75),aspect = "iso")

### xyplot - different colors for different cuts of the data
xyplot(lat ~ long, groups=Depth, data = quakes)


### using a theme to change color scheme in both plot and auto.key:

myTheme <- simpleTheme(col=c("green","red","blue"), lty=1:2)

densityplot(~ gcsescore,data=Chem97, groups=score, plot.points=FALSE, ref=TRUE, auto.key=list(columns=3), par.settings=myTheme)


myhist <- histogram(~ gcsescore | score,data=Chem97)
class(myhist)
summary(myhist)
plot(myhist)
myhist  ### also does the plot

## inside scripts, might need to do the plot(myhist) command explicitly - if called via source the plot doesn't always appear.

##### multiple plots per page - need to explicitly plot the trellis objects and print in a separate command, so I can use the split and more arguments
### Use split and more arguments to plot (or print is the same) - see  ?print.trellis

plot (histogram(~ gcsescore | factor(score),data=Chem97),split=c(1,1,1,2),more=TRUE)

plot(densityplot(~ gcsescore,data=Chem97, groups=score, plot.points=FALSE, ref=TRUE, auto.key=list(columns=3)),split=c(1,2,1,2),more=FALSE)

####split	 a vector of 4 integers, c(x,y,nx,ny) , that says to position the current plot at the x,y position in a regular array of nx by ny plots. (Note: this has origin at top left)
### more	 A logical specifying whether more plots will follow on this page.

##### might need to learn more about formula objects

########  bwplot = boxplots

bwplot( factor(score) ~ gcsescore | gender, data=Chem97, xlab="Average GCSE score")

bwplot( gcsescore^2.34 ~ gender | factor(score), data=Chem97, varwidth=TRUE, layout=c(6,1), ylab="Transformed GCSE score")

bwplot( Days ~ log(FSC.H), gvhd10, panel=panel.violin, box.ratio=3, xlab="log(Forward scatter)", ylab="Days past transplantation")

###### look at default settings etc
trellis.par.get()   ##  a list of lists
names(trellis.par.get())
# [1] grid.pars   fontsize background panel.background 
# [5] clip  add.line  add.text  plot.polygon     
# [9] box.dot  box.rectangle     box.umbrella      dot.line         
#[13] dot.symbol  plot.line         plot.symbol       reference.line   
#[17] strip.background  strip.shingle strip.border      superpose.line   
#[21] superpose.symbol  superpose.polygon regions           shade.colors     
#[25] axis.line  axis.text         axis.components   layout.heights   
#[29] layout.widths box.3d            par.xlab.text     par.ylab.text    
#[33] par.zlab.text  par.main.text     par.sub.text     
### to look at one of those lists
trellis.par.get()[["axis.line"]]
### and to look inside all of those lists
lapply(trellis.par.get(), names)

## to set option(s), for the current device:
trellis.par.set("axis.line"=list(col="transparent"))



lattice.getOption()   ##  a list of lists

