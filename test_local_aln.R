## align locally each sequence to each
# initialize distance matrix
DNA_dist <- matrix(nrow = length(seqs), ncol = length(seqs))
rownames(DNA_dist) <- names(seqs)
colnames(DNA_dist) <- names(seqs)

for (i in 1:length(seqs)) {
  for (j in 1:length(seqs)) {
    seq <- as.character(as.character(seqs[i]))
    seq2 <- as.character(as.character(seqs[j]))
    pa <- pairwiseAlignment(seq, seq2, type = "local")
    sd <- as.integer(stringDist(c(as.character(pattern(pa)), as.character(subject(pa)))))
    # if length of the alignment is very short, use stringDist without aligning
    if (width(pattern(pa)) > 65) {
      DNA_dist[i, j] <- sd
      DNA_dist[j, i] <- sd
    } else {
      unalignedSd <- as.integer(stringDist(c(seq, seq2)))
      DNA_dist[i, j] <- unalignedSd
      DNA_dist[j, i] <- unalignedSd
    }
  }
}

ccf0_256 <- mapply(function(x,y){ccf(x$dtcr,y, plot = FALSE)[0]$acf[1]}, x = cds_we_256, y = ribo_we_256)

mapply(function(x,y){FCN}, x = as.character(as.character(seqs)), y = as.character(as.character(seqs)))

width(pattern(pairwiseAlignment(x,y,type="local"))

as.integer(stringDist(c(as.character(pattern(pairwiseAlignment(x,y,type="local"))), as.character(subject(pairwiseAlignment(x,y,type="local"))))))

## first width
wid <- function(x,y){
  width(pattern(pairwiseAlignment(x,y,type="local")))
}
w <- outer(X=seqs,Y=seqs,FUN="wid")

## alignedSd, unalignedSD
alignedSd <- function(x,y) {
  as.integer(stringDist(c(as.character(pattern(pairwiseAlignment(x,y,type="local"))), as.character(subject(pairwiseAlignment(x,y,type="local"))))))
}

unalignedSd <- function(x,y) {
  as.integer(stringDist(c(x,y)))
}

checkWidth <- function(x) {
  if (x >= 70) alignedSd() else unalignedSd
}

y = function(x) if (x > 12) arrayInd(a, dim(o)) else NA

y = function(x) if (x %% 3) x+1 else NA
m <- matrix(1:10, nrow=2)
m[] <- vapply(m, y, numeric(1))
