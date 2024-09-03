Hi there! Welcome to my thesis GitHub repository. This repository contains the code I used to generate the plots and statistics presented in my research.

The FullCohortOutput files contain the basic functions applied to the sample, including filtering the cohort to obtain the 11,918 sample that will be analyzed. The files contain the function used to simulate metabolites, the three penalized regression functions, and the simulation of metabolites themselves. The bootstrap and Z-test results are contained within these files as well, which are generating from comparing the outcome values withheld in the testing test to the outcome values generated using the data in the training set. These two files can also be edited to subset the sample by ethnicity and then carry out the same analysis as was done on the full cohort.

The FullCohortOutput_AUC.R file creates tables that show the number of bootstrap tests with p-values less than 0.05 conducted on the entire filtered sample.

The FullCohortOutput_Correlation,Plots.R file creates tables that show the number of Z-tests tests with p-values less than 0.05 conducted on teh entire filtered sample, as well as produces plots of the AUC and Spearman correlation coefficients generated at a given training set size.

The DiminishingReturnsOfMetabolites.R file contains the simulation of 1 to 100 metabolites. It contains the calculation of the change in the AUC values (slope) for every 5 unit increase in the number of simulated metabolites and the calculation of the change in the Spearman correlation coeffcient (slope) for every 5 unit increase in the number of simulated metabolites. These slope values are calculated across all three penalized regression techniques and saved in dataframes. 

The DiminishingReturnsEvaluation.R file contains the calculations to determine where we no longer see an increase in the slope in both the AUC values and the Spearman corrleation coefficient values when adding more metabolites (point of diminishing returns).

If you have any questions or concerns, please send me a message. Cheers!
