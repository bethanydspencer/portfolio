#Experimentation Calculators

This file contains calculators for sample size and statistical significance for online tests using fixed horizon (frequentist) testing. 
Fixed horizon testing means that the sample size must be calculated before the test and significance checked once only after the test is complete, else the results are invalid.

There are calculators for two different types of metric: rates/probabilities (excel file) and  means (R scripts). 
Because rates and probabilities follow a binomial distribution, the sample size and significance is proportional to the rate/probability itself. Means can follow other distributions and in online experiments are often not normally distributed. 
For this reason the sample size calculator in R uses an empircal calculation of variance. 

Multiple comparisons
The excel sheet uses the Bonferroni correction to account for multiple variants, but still assumes only one metric will be tested per experiment. 
The R scripts currently do not have the Bonferroni correction and so should be used only for 1 metric and 1 variant. 
If you repeat this calculation for multiple metrics in the experiment, or multiple variants, the false positive rate will be much higher than the desired 5%.
The solution to this is to move to False Discovery Rate control, which will be built into the script in future.
