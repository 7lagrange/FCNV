
# prior emission mean and variance
emissionDistributions<-function(len=3000000000, binsize=1000, purity=0.8){
  length<-len/binsize
  numreads<-20*length
  tstate<-c(0,1,2,3,4,0,1,2,3,4,0,1,2,3,4)
  nstate<-c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3)
  p<-(tstate/2)/length
  q<-(nstate/2)/length
  mean<-log(purity*p/q+1-purity)
  variance<-(1-p)/(numreads*p)+(1-q)/(numreads*q)
  variance[p==0]<-(1-q[p==0])/(numreads*q[p==0])
  return(data.frame(tstate, nstate, mean, variance))
}

# initialize HMM model with states, transition matrix , emission matrix
initHMM = function(States, startProbs=NULL, transProbs=NULL,emissionProbs=NULL)
{
  nStates    = length(States)
  S          = rep(1/nStates,nStates)
  T          = 0.5*diag(nStates) + array(0.5/(nStates),c(nStates,nStates))
  names(S)   = States
  dimnames(T)= list(from=States,to=States)
  if(!is.null(startProbs)){S[]  = startProbs[]}
  if(!is.null(transProbs)){T[,] = transProbs[,]}
  return(list(States=States,startProbs=S,transProbs=T,
      emissionProbs=emissionProbs))
}


# calculate the probablity of observation under the normal distribution
emission_p = function(hmm, state,observation){
  mean_e <- hmm$emissionProbs[state,]$mean
  std_e <- sqrt(hmm$emissionProbs[state,]$variance)
  return (dnorm(observation,mean_e,std_e,log = TRUE)) }


# using viterbi algorithm to get the hidden states under parameter
viterbi = function(hmm, observation)
{ print (hmm)
  hmm$transProbs[is.na(hmm$transProbs)]  = 0 # remove NA in trans matrix
  nObservations  = length(observation)
  nStates    = length(hmm$States)
  v          = array(NA,c(nStates,nObservations))
  dimnames(v)= list(states=hmm$States,index=1:nObservations)
  # Init


  for(state in hmm$States)
  {
    v[state,1] = log(hmm$startProbs[state]) + emission_p(hmm, state,observation[1])
  }
  
  # Iteration
  for(k in 2:nObservations)
  {
    for(state in hmm$States)
    {
      maxi = NULL
      for(previousState in hmm$States)
      {
        temp = v[previousState,k-1] + log(hmm$transProbs[previousState,state]) 
        maxi = max(maxi, temp)
      }
      v[state,k] = emission_p(hmm, state,observation[k]) + maxi
    }
    if (all(is.na(v[,k])) ) {
    #if (v[1,k] == Inf ) {
    #if (all(v[,k] == rep(Inf,15)) ) {
      print (v[,k])
      print (v[,k-1])
      print (c(k,'NUMBER OF OBERSERVATION'))
      print (emission_p(hmm, state,observation[k]))
      break
    }
    #else{
    #  print (c(k,observation[k]))
    #  print (v[,k])
    #  break}
    
  }
  
  # Traceback
  viterbiPath = rep(NA,nObservations)
  print (v[,nObservations])
  for(state in hmm$States)
  {
    if(max(v[,nObservations])==v[state,nObservations])
    {
      viterbiPath[nObservations] = state
      break
    }
  }
  for(k in (nObservations-1):1)
  {
    for(state in hmm$States)
    {
      if(max(v[,k]+log(hmm$transProbs[,viterbiPath[k+1]]))
          ==v[state,k]+log(hmm$transProbs[state,viterbiPath[k+1]]))
      {
        viterbiPath[k] = state
        break
      }
    }
  }
  return(viterbiPath)
}



# using viterbi algorithm to get the most probable paramater
viterbiTraining = function(hmm, observation, maxIterations=100, delta=1E-9, pseudoCount=0)
{
  tempHmm = hmm
  tempHmm$transProbs[is.na(hmm$transProbs)]       = 0
  diff = c()
  for(i in 1:maxIterations)
  {
    # Counts
    vt = viterbiTrainingRecursion(tempHmm, observation)
    T  = vt$TransitionMatrix
    E  = vt$EmissionMatrix
    # Relativ Frequencies
    T[is.na(T)] = 0
    T = T/(apply(T,1,sum))
    
    d = sqrt(sum((tempHmm$transProbs-T)^2)) +  sqrt(sum((tempHmm$emissionProbs$mean-E$mean)^2))
    #print (d)
    diff = c(diff, d)
    tempHmm$transProbs    = T
    tempHmm$emissionProbs = E
    if(d < delta)
    {
      break
    }
  }
  tempHmm$transProbs[is.na(hmm$transProbs)]       = NA
  return(list(hmm=tempHmm,difference=diff))
}

# update the parameter 
viterbiTrainingRecursion = function(hmm, observation)
{
  TransitionMatrix    = hmm$transProbs
  TransitionMatrix[,] = 0
  EmissionMatrix      = hmm$emissionProbs

  v = viterbi(hmm,  observation)
  print (c(length(v),'asdasda'))
  for(i in 1:(length(observation)-1))
  { TransitionMatrix[v[i],v[i+1]] = TransitionMatrix[v[i],v[i+1]] + 1 }
  TransitionMatrix[TransitionMatrix==0]  = TransitionMatrix[TransitionMatrix==0] + 0.1#add psudo count
  
  
  for(i in  hmm$States){ # update emission matrix
  value <-   observation[v==i]
  mu = mean(value)
  sigma = var(value)

  
  n = length(value)
  if (n > 1){  
    if (sigma == 0) {sigma <- 0.1} # state without varance,set default
    m = EmissionMatrix$mean[i]
    tau = EmissionMatrix$variance[i]    
    EmissionMatrix$mean[i] =(n*tau*mu+sigma*m)/(n*tau+sigma)   
    EmissionMatrix$variance[i] = (2*sigma * tau)/(tau+sigma)  
  }
  

  }
  
  
  #print (TransitionMatrix)
  return(list(TransitionMatrix=TransitionMatrix,EmissionMatrix=EmissionMatrix))
}

