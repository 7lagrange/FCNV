**Clonal Fraction Hidden Markov Model**

We used hidden markov model to infer the tumor / normal cell admixure and CNV state.

All the programs are running on R. The code folder contains all the relevant scripts to generate the simulated data, process the data and post-process. 

**Functions**
functions in code/WGSsimulator.R are used to generate simulated data, you need to choose the simulated length, slide-bin size, the tumor purity options and tumor, as well as thr normal state transition probabilities (transition equallly likely to any different state). 
output will be in the data/ folder with filename ends with .simulated.csv

functions in code/find_states.R are used on the simulated data, it will apply the continuous HMM to the observed log.ratio and infer the hidden states at each position. You need to determine the maximum number of iterations. It takes several hours to iterate once on a regular PC with a full 3GB genome in 1kb bins.
output will be the in the data/ folder with two extra columns that indicate the hidden tumor states and normal states, filename will be end with .hmm.csv

functions in code/postprocess.R  are used to derive columns for posterior purity and accuracy of classification (with and without state compression) on the continuous HMM output.

**Example:**

open  example.R , execute all the command 

we run a small example on 30Mb data with 1 Kb bin size and two purity options of 0.9 and 0.8.

the output files are in the data folder:
threeK_90_80.simulated.csv contains simulated data from WGSsimulator.r
threeK_90_80.hmm.csv contains HMM output results from find_states.r
threeK_90_80.post.csv contains postprocess results from postprocess.r


