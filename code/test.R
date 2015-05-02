rm(list = ls())
source('HMM_modified.r')
simutation <- read.csv(file="simulation/threeM_90_90.csv",head=TRUE,sep=",")
simutation <- simutation[!is.na(simutation$Nstates),]
# Initial HMM
init_emission <- emissionDistributions()
init_transition <- matrix(rep(0.5/15,15*15),15,15) + 0.5*diag(15)
observation <- simutation$log.ratios
hmm = initHMM(c(1:15),
              transProbs=init_transition,
              emissionProbs=init_emission)



vt = viterbiTraining(hmm,observation,8)
print(vt)
hidden_state = viterbi(vt$hmm,observation)



summary(simutation$state-hidden_state)
hidden_Nstate <- hidden_state %/% 5 +1
hidden_Tstate <- (hidden_state-1) %% 5

simutation$hiddenstate <- hidden_state
simutation$hiddenNstate <- hidden_Nstate
simutation$hiddenTstate <- hidden_Tstate

correct <- simutation$state == hidden_state
simutation$consistance <- correct
write.csv(simutation, file = "simulation/threeM_90_90_output.csv")
