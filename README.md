# lability-critique
Citation to the publication associated with the dataset and code:

TO BE ADDED LATER

Author contact information: 

TO BE ADDED LATER

Brief summary of the published paper:

The goals of the paper were to: 
1.	Provide a clear, comprehensive protocol for demonstrating demographic lability (i.e., a higher stochastic population growth rate when vital rates vary over time than when they are fixed, due to non-linear relationships between vital rates and environmental drivers and/or skewed driver distributions)
2.	Use the protocol to assess whether any published studies show comprehensive evidence of lability
3.	Examine how lability might arise initially through natural selection acting on the shape of a vital rate function, relative to other life history adaptations to environmental stochasticity
4.	Examine whether selection would maintain lability once it has arisen
5.	Question previous predictions about the relationship between life history speed and the prevalence of demographic lability
6.	Examine types of environmental drivers that are – and are not – likely to support demographic lability
   
Description of the items in the repository:

Apart from a re-analysis of published data from Barraquand et al. 2014 (J. Animal Ecology 83: 375-387), the paper is theoretical, and this repository provides the code used to generate the figures in the paper and to perform a reanalysis of the data from Barraquand et al. 2014.   

Code was written by XXXXXXXXX, who also performed digitization and model fitting for the skua data from Barraquand et al. 2014.  Scripts were run in RStudio using R version 4.2.2 (2022-10-31 ucrt) running under Windows 10 x64 (build 19045).  As scripts were written for a PC platform, windows() commands are used to open graphical windows (they could be replaced with quartz() on a MAC machine).  Scripts are licensed under the MIT License described here: https://choosealicense.com/licenses/mit/

The script < fit skua productivity functions.R > uses data on skua productivity vs lemming data obtained by digitizing Fig. 2 in Barraquand et al. (2014). The digitized data are found in the file < skua data from Barraquand.csv >.

The data file < skua data from Barraquand.csv > has the following columns:
site: name of the site in Greenland, either Traill or Zackenberg
lemmings: average number of lemmings per hectare
fledglings: average number of fledglings per nest

The remaining R scripts in this repository create the figures in the main text and the supplement of the paper. The names of the scripts correspond to the figure numbers in the text and supplement of the manuscript. These scripts can be run in any order. 

The following is a list of these scripts and a brief description of the figures they make:

make fig 1.R - plots selection gradients for changes in the intercept, slope, and curvature of the adult survival function for a slow life history.

make fig 2.R – plots possible changes in a vital rate function under the persistent action of natural selection, potentially leading to the erosion of demographic lability

make fig 3.R – shows how adaptation to the local environment, through a unimodal vital rate function, could prevent demographic lability for fast and slow life histories

make fig S1-1.R – plots combinations of shape of a vital rate function and skewness in the distribution of the environmental driver that cause the vital rate’s arithmetic mean to be greater when the driver varies than when the driver is fixed at its median

make fig S1-2.R – plots the combinations of shape of a vital rate function and skewness in the distribution of the environmental driver at which the stochastic growth rate is higher when adult survival in the slow life history varies compared to when it is fixed at its value in either the mean or the median environment

make fig S4-1 S4-2.R – makes two plots showing local selection gradients acting on small portions of a vital rate function, one for a quadratic vital rate function and a symmetrical driver distribution (which produces lability at the driver mean), and one for a linear vital rate function and a skewed driver distribution (which produces lability relative to the driver median, but not the mean)

make fig S5-1.R – makes plots showing that the conclusion reached by LeCoeur et al. 2022 that demographic lability is more likely for species with a fast life history was predetermined by the assumptions they made in designing their simulations

make fig S6-1.R – plots the effect on the arithmetic mean vital rate and the stochastic growth rate of increasing mismatch between the optimal environment for a unimodal vital rate function and the mean of the (symmetrical) environmental distribution, for a range of values of the environmental standard deviation

make fig S7-1 S7-2.R – makes two plots of the effect of the shape of the adult survival function for a slow life history under a symmetric, bounded driver distribution on the stochastic growth rate, one plot for a quadratic vital rate function and one plot for adult survival following a power function  

make fig S7-3.R – effects of a skewed prey abundance distribution on the arithmetic mean of predator vital rates when predation follows a type III functional response (which is convex over  some driver range but concave over a different driver range)
