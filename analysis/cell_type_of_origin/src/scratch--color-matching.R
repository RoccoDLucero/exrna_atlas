

#Some test data
dat <- data.frame(x=runif(10),y=runif(10),
                  grp = rep(LETTERS[1:5],each = 2),stringsAsFactors = TRUE)

#Create a custom color scale
library(RColorBrewer)
myColors <- brewer.pal(5,"Set1")
names(myColors) <- levels(dat$grp)
colScale <- scale_colour_manual(name = "grp",values = myColors)

#One plot with all the data
p <- ggplot(dat,aes(x,y,colour = grp)) + geom_point()
p1 <- p + colScale

#A second plot with only four of the levels
p2 <- p %+% droplevels(subset(dat[4:10,])) + colScale




ggplot(subdata, aes(x = x, y = y, colour = fCategory)) +
    geom_point() +
    scale_colour_discrete(drop=TRUE,
                          limits = levels(dataset$fCategory))
