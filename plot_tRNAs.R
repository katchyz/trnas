## tRNA plots
library(plotrix)

d <- data.frame(a = c(1,2,3), b = c(4,5,6), c = c(7,8,9))
color2D.matplot(d)
# alignments (examples at least)
# colour-coded matrix

############
library(reshape)
library(ggplot2)
library(scales)

data <- structure(list(people = structure(c(2L, 3L, 1L, 4L), 
                                          .Label = c("bill", "mike", "sue", "ted"), 
                                          class = "factor"), 
                       apple = c(1L, 0L, 3L, 1L), 
                       orange = c(0L, 0L, 3L, 1L), 
                       peach = c(6L, 1L, 1L, 0L)), 
                  .Names = c("people", "apple", "orange", "peach"),
                  class = "data.frame", 
                  row.names = c(NA, -4L))
data.m <- melt(data)
data.m <- ddply(data.m, .(variable), transform, rescale = rescale(value))
p <- ggplot(data.m, aes(variable, people)) + 
  geom_tile(aes(fill = rescale), colour = "white") 
p + scale_fill_gradient(low = "white", high = "steelblue")

##############
distances <- read.csv("./trna_split/20aa/distances.csv")
LAwidth <- read.csv("./trna_split/20aa/LAwidth.csv")

distances.m <- melt(distances)
#ggplot(distances.m, aes(variable, X)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue")
ggplot(distances.m, aes(variable, X)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 8, breaks = c(1,8,80)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("distances")
ggsave(file = "/Users/kasia/Documents/PhD/scripts/tRNAs/plots/distances.png")

LAwidth.m <- melt(LAwidth)
#ggplot(LAwidth.m, aes(variable, X)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue")
ggplot(LAwidth.m, aes(variable, X)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 15, breaks = c(1,15,80)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("local alignment width")
ggsave(file = "/Users/kasia/Documents/PhD/scripts/tRNAs/plots/LAwidth.png")

## plot Gly_XXX 40-47
glyd <- distances[40:47,41:48]
glyd$X <- colnames(glyd)
glyd.m <- melt(glyd)

ggplot(glyd.m, aes(variable, X)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 8, breaks = c(1,8,80)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("distances, Gly")
ggsave(file = "./plots/gly_distances.png")

glyw <- LAwidth[40:47,41:48]
glyw$X <- colnames(glyw)
glyw.m <- melt(glyw)

ggplot(glyw.m, aes(variable, X)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 15, breaks = c(1,15,80)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("local alignment width, Gly")
ggsave(file = "./plots/gly_LAwidth.png")
