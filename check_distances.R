distances <- read.csv(file = "/Users/kasia/Documents/PhD/scripts/tRNAs/trna_split/20aa/distances.csv", header = TRUE)
distances$X <- NULL

LAwidth <- read.csv(file = "/Users/kasia/Documents/PhD/scripts/tRNAs/trna_split/20aa/LAwidth.csv", header = TRUE)
LAwidth$X <- NULL

p <- 0
for (i in 1:nrow(LAwidth)) {
  if (sum(LAwidth[i,][-i] > 14) > 1) {
    p <- p + 1
  }
}

for (i in 1:nrow(LAwidth)) {
  if (sum(LAwidth[i,][-i] > 15) <= 1) {
    print(i)
  }
}

# > names(LAwidth)[39]
# [1] "Glu_TTC2"