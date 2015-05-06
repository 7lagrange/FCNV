# FCNV

We used hidden markov model to infer the tumor / normal cell admixure and 

All the programs are running on R. The code code folder contain all the relevent script to generate the simulated data, process the data and evalutaion. 

Run code/WGSsimulator.R to get the simulated data, you need to choose the simulated length, slide-bin size as well as the tumor purity. 
output will be in the data/ folder.

Run code/find_states.R on the simulated data, it will apply the continuous HMM to the observed log.ratio and infer the hidden states at each position. You need to determine the maximum number of iteration. It takes hours to iterate once on a regular PC.
output will be the in the data/ folder with two extra column indicate the hidden tumor states and normal states.

Run code/postprocess.R to evaluate on the continuous HMM output, we mainly evaluate on the accuracy of the compressed data and purity.



