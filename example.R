rm(list = ls())

source('code/WGSsimulator.R')
source('code/HMM_modified.r')
source('code/find_states.R')
source('code/postprocess.R')


# parameter for simulate small size data of 3MB 
filename = 'threeK_90_80'
len=3*10^7
bin = 1000
tstatechange=0.0005
purity=c(0.9, 0.8)

# simulate data
WGSsimulator(filename,len=len,bin = bin,tstatechange=tstatechange, purity=purity)

# run HMM on data
find_state(filename)

# postprocess on data
postprocess(filename)
