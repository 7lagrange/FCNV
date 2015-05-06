# takes log.ration at each posiiton, using HMM to infer its  hidden Tumor and normal cell states 


rm(list = ls())
source('HMM_modified.r')

# read in data
simutation <- read.csv(file="data/threeM_90_90.csv",head=TRUE,sep=",")


# Initial HMM
init_emission <- emissionDistributions()
init_transition <- matrix(rep(0.5/15,15*15),15,15) + 0.5*diag(15)
observation <- simutation$log.ratios
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
write.csv(simutation, file = "data/threeM_90_90_output.csv")
