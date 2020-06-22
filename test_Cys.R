
seqs <- readDNAStringSet('Cys_nospaces.txt')
clusterSequences(seqs,max.Clusters=5)

##Run functions
clusterSequences = function(seqs, tryClusterRange = 20, maxEdits = 8, threshold=0.75, max.Clusters=NULL){
  
  # define recursive function
  clusterSeqsRecursive = function(seqs,tryClusterRange,maxEdits){
    
    library(Biostrings)
    library(DECIPHER)
    
    
    if(length(seqs)==1){return(seqs)}
    
    getWss <- function(d) {
      sum(scale(d, scale = FALSE)^2)
    }
    
    getWithinClusterSum <- function(hc, k, d) {
      d<- as.matrix(d)
      cl <- cutree(hc, k)
      spl <- split(d,cl)
      wss <- sum(sapply(spl, getWss))
      wss
    }
    
    # alternative option: use kmeans and increase the nstart parameter 
    DNA <- AlignSeqs(seqs)
    DNA_dist<- as.dist(as.matrix(stringDist(seqs)))
    DNA_hclust<- hclust(DNA_dist)
    
    if(!is.null(max.Clusters)) {
      # force max numbers of clusters
      message(paste("forcing number of clusters to be:",max.Clusters))
      warning("ignoring parameters threshold and maxEdits")
      clIds<- cutree(DNA_hclust,k=max.Clusters)
      clNam<- names(clIds)
      clusters<-list()
      for (cl in unique(clIds)) {
        clusters[[cl]]<- seqs[clNam[clIds==cl]]
      }
      ###
      consensusList<- list()
      i <- 1
      for(cl in clusters){
        clIdsX <- names(cl)
        seqsX <- seqs[clIdsX]
        consensusList[i] <- consensusString(seqsX, ambiguityMap="N",threshold=threshold)
        i <- i + 1
      }
      ###
      tooManyMismatches <- list()
      i <- 1
      for(cl in clusters){
        cons<- DNAStringSet(consensusList[[i]])
        input<- c(as.character(as.character(cl)),as.character(cons))
        res<-as.matrix(stringDist(input))
        distances<-res[nrow(res),]
        print(consensusList[[i]])
        print(table(distances))
        tooManyMismatches[[i]] <- names(input)[distances > 8] #contains list with outliers
        i <- i + 1
      }
      ###
    } else {
      # recursive clustering
      res<- c()
      krange<- min(tryClusterRange,length(seqs))
      for(k in 1:krange){
        res<- c(res, getWithinClusterSum(DNA_hclust, k, DNA_dist))
      }
      #barplot(res, names.arg = 1:krange)
      
      #kmeans
      
      # get optimal cluster size
      for(i in seq(2,length(res), 1)){ if(res[i]>=0.99*res[i-1]){break} }
      
      print(i)
      clIds<- cutree(DNA_hclust, i)
      clNam<- names(clIds)
      
      # make a consensus list
      consensusList<- list()
      for(cl in clIds){
        clIdsX<- clNam[clIds==cl]
        seqsX<- seqs[clIdsX]
        consensusList[[cl]]<- consensusString(seqsX, ambiguityMap="N",threshold=threshold)
      }
      
      # recursion
      clusters<-list()
      for (cl in unique(clIds)) {
        # if the distances are > maxEdits, call 
        input<- c(seqs[clNam[clIds==cl]],DNAStringSet(consensusList[[cl]]))
        res<- as.matrix(stringDist(input))
        distances<- res[nrow(res),-ncol(res)]
        inputNames<- names(input)[-ncol(res)]
        print(consensusList[[cl]])
        print(table(distances))
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
  message(paste("found",length(res),"clusters"))
  message(paste("cluster sizes:"))
  message(paste(sapply(res,length),collapse=","))
  return(res)
  
}