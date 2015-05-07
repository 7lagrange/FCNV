#uses markov chain (efficiently simulated by geometric variables, given state change probabilities for tumor and normal)
#Nstates 1,2,3
#Tstates 0,1,2,3,4



readStarts<-function(length, numreads, CNV.data){
  CNV.data$Ncount<-(tabulate(sample(1:length, numreads, replace=TRUE, prob=CNV.data$Nstates)))
  CNV.data$Tcount<-(tabulate(sample(1:length, numreads, replace=TRUE, prob=CNV.data$Tstates)))  
  CNV.data$log.ratios<-log((CNV.data$tpurity*CNV.data$Tcount+(1-CNV.data$tpurity)*CNV.data$Ncount)/CNV.data$Ncount)      
  return(CNV.data)                               
}

SimulateData<-function(len=300000000, binsize=10000, tstatechange=0.0005, nstatechange=0.0001, purity=c(0.9, 0.7)){
  length<-len/binsize
  numreads<-20*length
  Nstates<-rep(2, length)
  Tstates<-rep(2, length)
  tpurity<- rep (purity[1], length)
  Ncount<- rep(0, length)
  Tcount<- rep(0, length)
  log.ratios<- rep(0, length)
  CNV.data<-data.frame(Nstates, Tstates, tpurity, Ncount, Tcount, log.ratios)
    
  tindex<-rgeom(1, tstatechange)
  tnext<-tindex+rgeom(1, tstatechange)
  nindex<-rgeom(1, nstatechange)
  nnext<-nindex+rgeom(1,nstatechange)
  while(tnext<length){
    state<-sample(setdiff(0:4, CNV.data$Tstates[tindex-1]),1)
    CNV.data$Tstates[tindex:tnext]<-state
    CNV.data$tpurity[tindex:tnext]<-sample(purity,1)
    tindex<-tnext
    tnext<-tnext+rgeom(1,tstatechange)
  }
  while(nnext<length){
    state<-sample(setdiff(1:3, CNV.data$Nstates[nindex-1]),1)
    CNV.data$Nstates[nindex:nnext]<-state
    nindex<-nnext
    nnext<-nnext+rgeom(1,nstatechange)
  }
  
  CNV.data<-readStarts(length, numreads, CNV.data)
  return(CNV.data)
}


WGSsimulator<-function(output,len=3*10^9,bin = 1000,tstatechange=0.0005, purity=c(0.9, 0.8)){
  threeK <- SimulateData(len,bin ,tstatechange, purity)
  
  # Tumor/Normal state transformation 
  threeK$state <- threeK$Tstates + 1 + (threeK$Nstates-1)*5

  # remove invalid data
  threeK <- threeK[threeK$log.ratios!= Inf,]
  threeK <- threeK[!is.na(threeK$Nstates),]
  
  
  # write to csv for HMM and post-process
  output <- paste("data/",output,".simulated.csv",sep="")
  write.csv(threeK, file =output ,row.names=FALSE)  
  
}




