# A script in R to perform an Random Forest Approximate Bayesian Computation analysis under the Isolation with Migration model from microsatellite data
Script in R to perform an Approximate Bayesian Computation (ABC) analysis under the Isolation with Migration model from microsatellite data.

Demographic history was inferred from microsatellite data using an approximate Bayesian computation (ABC) approach via random forests (Marin et al. 2016; Pudlo et al. 2016). In this approach, data is simulated from the demographic model with parameter values taken from prior probability distributions and data is transformed into summary statistics. Random forests are used to learn about the parameters from the simulated summary statistics. The resulting random forests can then be used to estimate the posterior probability distributions of parameters.

In this analysis, a model of two populations (eastern and western clusters) is evaluated. Each population is characterized by a parameter θ (θW=4NWμ and θE=4NEμ, where NW is the effective population size of the western population, NE is the effective population size of the eastern population, and μ is the mutation rate). Western population is founded by individuals from the eastern population at time T=t/4NW (time tF measured in number of generations). Two concurrent models are evaluated regarding the presence or absence of gene flow between the two populations. In the case of presence of gene flow there is an additional parameter, the scaled migration rate M=4NWm. Microsatellites are assumed to mutate following a generalized stepwise mutation model (GSM), in which the number of repeat units gained or lost in each mutation is taken from a geometric distribution with parameter PGSM. Data under this model is generated by simulation using coalescent simulator ms (Hudson 2002) with a custom script (see below) to transform its output into microsatellite data. Each simulated data set was summarized by statistics used in population genetics to characterize microsatellite genetic diversity and population differentiation, known to be informative about demography.

Data input: file in STRUCTURE format (.str) with alleles coded in number of repeats (microsatellite data), two rows per individual, second column indicating the population, missing data coded as -9

### REFERENCES

Hudson, R. R. 2002. Generating samples under a Wright-Fisher neutral model of genetic variation. Bioinformatics 18:337–338.

Marin, J.-M., L. Raynal, P. Pudlo, M. Ribatet, and C. P. Robert. 2016. ABC random forests for Bayesian parameter inference. arXiv:1605.05537.

Pudlo, P., J.-M. Marin, A. Estoup, J.-M. Cornuet, M. Gautier, and C. P. Robert. 2016. Reliable ABC model choice via random forests. *Bioinformatics* **32**:859–866.


### Requirements

Executable of ms and script to read ms output (distributed woth ms) in /bin

R packages: plyr, pegas, mmod, hexbin, grid, abcrf, quantregForest


