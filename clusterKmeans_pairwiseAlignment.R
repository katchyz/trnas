library(DECIPHER)
library(Biostrings)

path <- "./trna_split/20aa"
#path <- "/export/valenfs/projects/Zoya/tRNA_microarray/tRNA/split/20aa"
files <- list.files(path=path, pattern="*.txt", full.names=T, recursive=FALSE)

l <- lapply(files, function(x) {
  seqs <- readDNAStringSet(x)
  res <- clusterSequencesKmeans(seqs) ## ??
  summarize_output(res)
})

outputs <- list()
ress <- list()
for (file in files) {
  seqs <- readDNAStringSet(file)
  res <- clusterSequencesKmeans(seqs)
  #summarize_output(res)
  aa <- substr(file, nchar(file)-10, nchar(file)-4)
  ress[[aa]] <- res
  outputs[[aa]] <- summarize_output(res)
}

save(ress, file = "./trna_split/20aa/ress.Rsave")
save(outputs, file = "./trna_split/20aa/outputs.Rsave")

###Functions to run
clusterSequencesKmeans = function(seqs, tryClusterRange = 20, maxEdits = 8, threshold=0.5, max.Clusters=NULL){
  
  # define recursive function
  clusterSeqsRecursive = function(seqs,tryClusterRange,maxEdits){
    
    
    if(length(seqs)==0){return(NULL)}
    if(length(seqs)==1){return(seqs)}
    
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
    
    if(!is.null(max.Clusters)) {
      # force max numbers of clusters
      message(paste("forcing number of clusters to be:",max.Clusters))
      warning("ignoring parameters threshold and maxEdits")
      DNA_km<- kmeans(DNA_dist,centers=max.Clusters,nstart=50)
      clIds<- DNA_km$cluster
      clNam<- names(clIds)
      clusters<-list()
      for (cl in unique(clIds)) {
        clusters[[cl]]<- seqs[clNam[clIds==cl]]
      }
    } else {
      # recursive clustering
      res<- c()
      kmClusters<- list()
      krange<- min(tryClusterRange,length(seqs)-1)
      print(length(seqs))
      #print(krange)
      for(k in 1:krange){
        km<- kmeans(DNA_dist,centers=k,nstart=50)
        kmClusters[[k]]<- km
        res<- c(res, km$betweenss/km$totss)
      }
      
      # get optimal cluster size
      if(length(res)==1){
        return(seqs)
      } else {
        message("find optimal number of clusters...")
        for(i in seq(2,length(res), 1)){ if((res[i]-res[i-1])/res[i-1] <= 0.1){break} }
      }
      
      print(i)
      clIds<- kmClusters[[i]]$cluster
      print(table(clIds))
      clNam<- names(clIds)
      
      # make a consensus list
      consensusList<- list()
      for(cl in clIds){
        clIdsX<- clNam[clIds==cl]
        seqsX<- seqs[clIdsX] ### ?????? should it be on pairwise aligned ??????????
        consensusList[[cl]]<- consensusString(seqsX, ambiguityMap="N",threshold=threshold)
      }
      
      # recursion
      clusters<-list()
      for (cl in unique(clIds)) {
        # if the distances are > maxEdits, call 
        input<- c(seqs[clIds==cl],DNAStringSet(consensusList[[cl]]))
        res<- as.matrix(stringDist(input))
        distances<- res[nrow(res),-ncol(res)]
        inputNames<- names(input)[-ncol(res)]
        
        tooManyMismatches<- inputNames[distances > maxEdits]
        # call the function again on a subset (recursion)
        if(length(tooManyMismatches)>0){
          message("recluster...")
          clusters[[cl]]<- clusterSeqsRecursive(seqs[clNam[clIds==cl]],tryClusterRange,maxEdits)   
        }else {
          clusters[[cl]]<- seqs[clNam[clIds==cl]]
        }
      }
    }
    
    return(clusters)
  }
  
  # run
  res<-clusterSeqsRecursive(seqs=seqs,tryClusterRange = tryClusterRange,maxEdits = maxEdits)
  res<- unlist(res,recursive = TRUE)
  res<- res[order(sapply(res,length),decreasing=TRUE)]
  message(paste("found",length(res),"clusters"))
  message(paste("cluster sizes:"))
  message(paste(sapply(res,length),collapse=","))
  return(res)
}
summarize_output <- function(res){
  consensusList<- list()
  for(i in 1:length(res)){
    #consensusList[[i]] <- consensusString(consensusMatrix(res[[i]]),ambiguityMap="N")
    consensusList[[i]]<- ConsensusSequence(res[[i]],threshold=0.5) ### ?????? should it be on pairwise aligned?
    ### !!!!!!!!!! fix ambiguity mapping
  }
  tooManyMismatches <- list()
  for(i in 1:length(res)){
    cons<- DNAStringSet(consensusList[[i]])
    input<- c(res[[i]],cons)
    res2<-as.matrix(stringDist(input))
    distances<-res2[nrow(res2),]
    print(consensusList[[i]])
    print(table(distances))
    tooManyMismatches[[i]] <- names(input)[distances > 8] #contains list with outliers
  }
  
  #total number of mismatches
  tooManymismatches2 <- list()
  tooManyMismatches<- list()
  for(i in 1:length(res)){
    cons<- DNAStringSet(consensusList[[i]])
    for(j in 1:length(res)){
      input<- c(res[[j]],cons)
      res2<-as.matrix(stringDist(input))
      distances<-res2[nrow(res2),]
      #print(consensusList[[i]])
      #print(table(distances))
      tooManyMismatches[[j]] <- names(input)[distances > 8] #contains list with outliers
    }
    tooManymismatches2[[i]] <- tooManyMismatches 
  }
  
  lapply(tooManymismatches2[[1]],length) #mismatches to each consensus seq per cluster (for cluster 1, 12 mismathces)
  sum_mismatch <- list()
  
  for(i in 1:length(res)){
    sum_mismatch[[i]] <- lapply(tooManymismatches2[[i]],length)
  }
  mismatch_mat <- matrix(unlist(sum_mismatch), ncol = length(res), byrow = TRUE)
  mismatch_min <- c()
  for(i in 1:length(res)){
    mismatch_min[i]<- min(mismatch_mat[,i])
  }
  total_mismatches <- sum(mismatch_min)
  message(paste("Total number of mismatches:",total_mismatches))
  message("Total number of sequences analyzed: ",length(seqs))
  return(consensusList)
}


