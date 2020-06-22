# cluster together for each anticodon ?

for (anticodon in outputs) {
  if (length(anticodon) == 2) {
    print(pairwiseAlignment(anticodon[[1]], anticodon[[2]], type = "local"))
  }
}

# check all against all - if they would hybridize together (pairwiseAlignment) ?
# make a matrix
u <- unlist(outputs)
uc <- sapply(u, function(x){as.character(as.character(x))})

LAwidth <- matrix(nrow=length(uc), ncol=length(uc))
rownames(LAwidth) <- names(uc)
colnames(LAwidth) <- names(uc)
for (i in names(uc)) {
  for (j in names(uc)) {
    pa <- pairwiseAlignment(uc[i],uc[j],type="local")
    LAwidth[i,j] <- width(pattern(pa))
    LAwidth[j,i] <- width(pattern(pa))
  }
}

m <- as.matrix(stringDist(uc))
write.csv(m, "/Users/kasia/Documents/PhD/scripts/tRNAs/trna_split/20aa/distances.csv")

write.csv(LAwidth, "/Users/kasia/Documents/PhD/scripts/tRNAs/trna_split/20aa/LAwidth.csv")

p <- 0
for (i in 1:nrow(m)) {
  if (sum(m[i,] < 8) > 1) {
    p <- p + 1
  }
}

# save consensus sequences