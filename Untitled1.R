##
library(Biostrings)

tRNAdir <- '/Volumes/USELESS/DATA/fasta/danio_rerio/trna_split'

seqs <- readDNAStringSet('/Volumes/USELESS/DATA/fasta/danio_rerio/trna_split/Ala_AGC.txt')

# library(seqinr)
# read.fasta()

n = 4
res <- clusterSequencesKmeans(seqs,max.Clusters=n)
output_Asn_ATT <- summarize_output(res) #not matching 2/15

tRNA <- "Asn"
anticodon <- "GTT"
seqs <- readDNAStringSet(paste(tRNA,"_",anticodon,".txt",sep=""))
n=2
#res_Asn_GTT <-clusterSequencesKmeans(seqs,max.Clusters=n)
res <- clusterSequencesKmeans(seqs,max.Clusters=n)
output_Asn_GTT <- summarize_output(res) #not matching 21/470

#matching function for Asn
tooManymatches <- list()
tooManymatches_seqs <- list()
tooManymatches2 <- list()
tooManymatches_seqs2 <- list()
ConsensusSequence(res[[i]], threshold=0.5,minInformation=0.5)
for(i in 1:4){ #n needs to be number of ATT seqs
  cons<- DNAStringSet(output_Asn_ATT[[i]]) #consensus list from Asn_ATT
  for(j in 1:n){
    input<- c(res[[j]],cons)
    res2<-as.matrix(stringDist(input))
    distances<-res2[nrow(res2),]
    print(output_Asn_ATT[[i]])
    print(table(distances))
    tooManymatches_seqs[[j]] <- input[distances < 8]
    tooManymatches[[j]] <- names(input[distances < 8]) #contains list with matching seq names
  }
  tooManymatches2[[i]] <- tooManymatches
  tooManymatches_seqs2[[i]] <- tooManymatches_seqs
}

#total matching sequences to other consensus seqs
sum_match <- list()
for(i in 1:4){ #change to number of ATT seqs
  sum_match[[i]] <- lapply(tooManymatches2[[i]],length)#matches to each consensus seq per cluster
}

match_mat <- matrix(unlist(sum_match), ncol = 4, byrow = TRUE) #change num
match_max <- c()
for(i in 1:4){ #change num
  match_max[i]<- sum(match_mat[,i])
}
total_matches <- sum(match_max)
print(total_matches) ##454

#matching of GTT consensus seqs to ATT clusters
tooManymatches <- list()
tooManymatches_seqs <- list()
tooManymatches2 <- list()
tooManymatches_seqs2 <- list()

for(i in 1:2){ #n needs to be number of GTT seqs
  cons<- DNAStringSet(output_Asn_GTT[[i]]) #consensus list from Asn_GTT
  for(j in 1:n){
    input<- c(res[[j]],cons)
    res2<-as.matrix(stringDist(input))
    distances<-res2[nrow(res2),]
    print(output_Asn_ATT[[i]])
    print(table(distances))
    tooManymatches_seqs[[j]] <- input[distances < 8]
    tooManymatches[[j]] <- names(input[distances < 8]) #contains list with matching seq names
  }
  tooManymatches2[[i]] <- tooManymatches
  tooManymatches_seqs2[[i]] <- tooManymatches_seqs
}

#total matching sequences to other consensus seqs
sum_match <- list()
for(i in 1:2){ #change to number of ATT seqs
  sum_match[[i]] <- lapply(tooManymatches2[[i]],length)#matches to each consensus seq per cluster
}

match_mat <- matrix(unlist(sum_match), ncol = 2, byrow = TRUE) #change num
match_max <- c()
for(i in 1:2){ #change num
  match_max[i]<- sum(match_mat[,i])
}
total_matches <- sum(match_max)
print(total_matches)

###ARG
tRNA <- "Arg"
anticodon <- "GCG"
seqs <- readDNAStringSet(paste(tRNA,"_",anticodon,".txt",sep=""))
n=2
res <- clusterSequencesKmeans(seqs,max.Clusters=n)
output_Arg_GCG <- summarize_output(res) #not matching 0/7

tRNA <- "Arg"
anticodon <- "TCG"
seqs <- readDNAStringSet(paste(tRNA,"_",anticodon,".txt",sep=""))
n=3
res <- clusterSequencesKmeans(seqs,max.Clusters=n)
output_Arg_TCG <- summarize_output(res) #not matching 3/50

tRNA <- "Arg"
anticodon <- "TCT"
seqs <- readDNAStringSet(paste(tRNA,"_",anticodon,".txt",sep=""))
n=2
res <- clusterSequencesKmeans(seqs,max.Clusters=n)
output_Arg_TCT <- summarize_output(res) #not matching 6/139


#matching function for Asn
tooManymatches <- list()
tooManymatches_seqs <- list()
tooManymatches2 <- list()
tooManymatches_seqs2 <- list()

for(i in 1:4){ #n needs to be number of ATT seqs
  cons<- DNAStringSet(output_Asn_ATT[[i]]) #consensus list from Asn_ATT
  for(j in 1:n){
    input<- c(res[[j]],cons)
    res2<-as.matrix(stringDist(input))
    distances<-res2[nrow(res2),]
    print(output_Asn_ATT[[i]])
    print(table(distances))
    tooManymatches_seqs[[j]] <- input[distances < 8]
    tooManymatches[[j]] <- names(input[distances < 8]) #contains list with matching seq names
  }
  tooManymatches2[[i]] <- tooManymatches
  tooManymatches_seqs2[[i]] <- tooManymatches_seqs
}

#total matching sequences to other consensus seqs
sum_match <- list()
for(i in 1:4){ #change to number of ATT seqs
  sum_match[[i]] <- lapply(tooManymatches2[[i]],length)#matches to each consensus seq per cluster
}

match_mat <- matrix(unlist(sum_match), ncol = 4, byrow = TRUE) #change num
match_max <- c()
for(i in 1:4){ #change num
  match_max[i]<- sum(match_mat[,i])
}
total_matches <- sum(match_max)
print(total_matches) ##454




#matching function
tooManymatches <- list()
tooManymatches_seqs <- list()
tooManymatches2 <- list()
tooManymatches_seqs2 <- list()
consensusList<- list()
for(i in 1:n){
  consensusList[[i]] <- consensusString(consensusMatrix(res[[i]]),ambiguityMap="N")
} 
for(i in 1:n){
  cons<- DNAStringSet(consensusList[[i]])
  #print(i)
  tooManymatches <- list()
  tooManymatches_seqs <- list()
  for(j in 1:n){
    if(i!=j){
      input<- c(res[[j]],cons)
      res2<-as.matrix(stringDist(input))
      distances<-res2[nrow(res2),]
      #print(consensusList[[i]])
      #print(table(distances))
      tooManymatches_seqs[[j]] <- input[distances < 8]
      tooManymatches[[j]] <- names(input[distances < 8]) #contains list with matching seq names
      #print(j)
    }
  }
  tooManymatches2[[i]] <- tooManymatches
  tooManymatches_seqs2[[i]] <- tooManymatches_seqs
}

#total matching sequences to other consensus seqs
sum_match <- list()
for(i in 1:n){
  sum_match[[i]] <- lapply(tooManymatches2[[i]],length)#matches to each consensus seq per cluster
}

match_mat <- matrix(unlist(sum_match), ncol = n, byrow = TRUE)
match_max <- c()
for(i in 1:n){
  match_max[i]<- sum(match_mat[,i])
}
total_matches <- sum(match_max)
print(total_matches -n*(n-1)) ##make sure this number is correct

#export list that match to other consensus sequences
tRNA = "Ala"
matches_seqs.df <- data.frame()
for(i in 1:n){
  matches_seqs.df <- do.call("rbind", lapply(tooManymatches_seqs2[[i]], as.data.frame)) 
  write.csv(matches_seqs.df, file = paste('match_seqs', tRNA,i,'csv', sep = '.')) 
}



###Functions to run
clusterSequencesKmeans = function(seqs, tryClusterRange = 20, maxEdits = 8, threshold=0.5, max.Clusters=NULL){
  
  require(Biostrings)
  require(DECIPHER)
  
  # define recursive function
  clusterSeqsRecursive = function(seqs,tryClusterRange,maxEdits){
    
    
    if(length(seqs)==0){return(NULL)}
    if(length(seqs)==1){return(seqs)}
    
    # alternative option: use kmeans and increase the nstart parameter 
    DNA_dist<- as.matrix(stringDist(seqs))
    #DNA_hclust<- hclust(DNA_dist)
    
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
        seqsX<- seqs[clIdsX]
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
  for(i in 1:n){
    #consensusList[[i]] <- consensusString(consensusMatrix(res[[i]]),ambiguityMap="N")
    consensusList[[i]]<- ConsensusSequence(res[[i]],threshold=0.5)
  }
  tooManyMismatches <- list()
  for(i in 1:n){
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
  for(i in 1:n){
    cons<- DNAStringSet(consensusList[[i]])
    for(j in 1:n){
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
  
  for(i in 1:n){
    sum_mismatch[[i]] <- lapply(tooManymismatches2[[i]],length)
  }
  mismatch_mat <- matrix(unlist(sum_mismatch), ncol = n, byrow = TRUE)
  mismatch_min <- c()
  for(i in 1:n){
    mismatch_min[i]<- min(mismatch_mat[,i])
  }
  total_mismatches <- sum(mismatch_min)
  message(paste("Total number of mismatches:",total_mismatches))
  message("Total number of sequences analyzed: ",length(seqs))
  return(consensusList)
}


