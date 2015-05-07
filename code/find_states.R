# takes log.ration at each posiiton, using HMM to infer its  hidden Tumor and normal cell states 
# read in data

find_state<-function(filename){

input <- paste("data/",filename,".simulated.csv",sep="")
simutation <- read.csv(file=input,head=TRUE,sep=",")
observation <- simutation$log.ratios


# Initial HMM
p <- mean(purityPrior(observation),na.rm=TRUE)
init_emission <- emissionDistributions(purity = p)
init_transition <- matrix(rep(0.5/15,15*15),15,15) + 0.5*diag(15)

hmm = initHMM(c(1:15),
              transProbs=init_transition,
              emissionProbs=init_emission)

# run continuous HMM
vt = viterbiTraining(hmm,observation,10)
print(vt)
hidden_state = viterbi(vt$hmm,observation)
summary(simutation$state-hidden_state)

# decode tumor/normal state from HMM output
hidden_Nstate <- hidden_state %/% 5 +1
hidden_Tstate <- (hidden_state-1) %% 5
simutation$hiddenstate <- hidden_state
simutation$hiddenNstate <- hidden_Nstate
simutation$hiddenTstate <- hidden_Tstate

# compared hidden state to real state 
correct <- simutation$state == hidden_state
simutation$consistance <- correct

output <- paste("data/",filename,".hmm.csv",sep="")
write.csv(simutation, file = output)
}
