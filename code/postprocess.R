# postprocess on HMM output, calculate the purity and accuracy for compressed hidden states


# function to get postpurity
postp<-function(CNV.data,numclusters=2){
  CNV.data$postpurity<-(exp(CNV.data$log.ratio)-1)/(CNV.data$hiddenTstate/CNV.data$hiddenNstate-1)
  CNV.data$postpurity[CNV.data$hiddenTstate==CNV.data$hiddenNstate]<-NA
  CNV.data$postpurity[CNV.data$postpurity>1]<-NA
  CNV.data$postpurity[CNV.data$postpurity<0]<-NA
  clust<-kmeans(CNV.data$postpurity[!is.na(CNV.data$postpurity)], centers=numclusters)
  CNV.data$clusterpurity<-CNV.data$postpurity
  CNV.data$clusterpurity[!is.na(CNV.data$postpurity)]<-clust$centers[clust$cluster]
  CNV.data$truepuritydiff<-CNV.data$clusterpurity
  CNV.data$truepuritydiff[!is.na(CNV.data$clusterpurity)]<-CNV.data$tpurity[!is.na(CNV.data$clusterpurity)]-CNV.data$clusterpurity[!is.na(CNV.data$clusterpurity)]
  return(CNV.data)
}

# function to get compress states 
compressStates<-function(CNV.data){
  CNV.data$compressedState<-CNV.data$state
  CNV.data$compressedState[CNV.data$state==1 | CNV.data$state==6 | CNV.data$state==11]<-1
  CNV.data$compressedState[CNV.data$state==7 | CNV.data$state==12 | CNV.data$state==13]<-2
  CNV.data$compressedState[CNV.data$state==2 | CNV.data$state==8 | CNV.data$state==14]<-3
  CNV.data$compressedState[CNV.data$state==9 | CNV.data$state==15]<-4
  CNV.data$compressedState[CNV.data$state==3 | CNV.data$state==10]<-5
  CNV.data$compressedState[CNV.data$state==4 | CNV.data$state==5]<-6
  
  CNV.data$compressedhiddenState<-CNV.data$hiddenstate
  CNV.data$compressedhiddenState[CNV.data$hiddenstate==1 | CNV.data$hiddenstate==6 | CNV.data$hiddenstate==11]<-1
  CNV.data$compressedhiddenState[CNV.data$hiddenstate==7 | CNV.data$hiddenstate==12 | CNV.data$hiddenstate==13]<-2
  CNV.data$compressedhiddenState[CNV.data$hiddenstate==2 | CNV.data$hiddenstate==8 | CNV.data$hiddenstate==14]<-3
  CNV.data$compressedhiddenState[CNV.data$hiddenstate==9 | CNV.data$hiddenstate==15]<-4
  CNV.data$compressedhiddenState[CNV.data$hiddenstate==3 | CNV.data$hiddenstate==10]<-5
  CNV.data$compressedhiddenState[CNV.data$hiddenstate==4 | CNV.data$hiddenstate==5]<-6

  CNV.data$compressedconsistence<-CNV.data$compressedState==CNV.data$compressedhiddenState
  return(CNV.data)
}

# function to cluster states using naive threshold
naivethreshold<-function(CNV.data){
  purity<-mean(purityPrior(CNV.data$log.ratios),na.rm=TRUE)
  normal<-log((purity*2+(1-purity)*2)/2)
  partdel<-log((purity*2+(1-purity)*3)/3)
  partamp<-log((purity*4+(1-purity)*3)/3)
  homoamp<-log((purity*4+(1-purity)*2)/2)
  homodel<-log((purity*0+(1-purity)*2)/2)
  largeamp<-log((purity*3+(1-purity)*1)/1)
  CNV.data$tcompressedstate<-CNV.data$compressedstate
  CNV.data$tcompressedstate[CNV.data$log.ratio<=homodel]<-1
  CNV.data$tcompressedstate[CNV.data$log.ratio>homodel & CNV.data$log.ratio<=partdel]<-2
  CNV.data$tcompressedstate[CNV.data$log.ratio>partdel & CNV.data$log.ratio<=partamp]<-3
  CNV.data$tcompressedstate[CNV.data$log.ratio>partamp & CNV.data$log.ratio<=homoamp]<-4
  CNV.data$tcompressedstate[CNV.data$log.ratio>homoamp & CNV.data$log.ratio<=largeamp]<-5
  CNV.data$tcompressedstate[CNV.data$log.ratio>=largeamp]<-6

  return(CNV.data)
}


postprocess <- function(filename){

input <- paste("data/",filename,".hmm.csv",sep="")

simutation <- read.csv(file=input,head=TRUE,sep=",")

# post process HMM output
result <- naivethreshold(postp(compressStates(simutation)))


output <- paste("data/",filename,".post.csv",sep="")
write.csv(result, file = output)

}

