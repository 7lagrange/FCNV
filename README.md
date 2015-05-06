# FCNV

We used hidden markov model to infer the tumor / normal cell admixure and CNV state 

All the programs are running on R. The code folder contains all the relevant scripts to generate the simulated data, process the data and post-process. 

Run code/WGSsimulator.R to get the simulated data, you need to choose the simulated length, slide-bin size as well as the tumor purity options and tumor, normal state transition probabilities (transition equallly likely to any different state). 
output will be in the data/ folder.

Run code/find_states.R on the simulated data, it will apply the continuous HMM to the observed log.ratio and infer the hidden states at each position. You need to determine the maximum number of iterations. It takes several hours to iterate once on a regular PC with a full 3GB genome in 1kb bins.
output will be the in the data/ folder with two extra columns that indicate the hidden tumor states and normal states.

Run code/postprocess.R to derive columns for posterior purity and accuracy of classification (with and without state compression) on the continuous HMM output.



