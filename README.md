# Ancestor regression
R Software for performing ancestor regression in structural vector autoregressive models

lin-anc-ts.R: Contains the main function for ancestor regression <br>
Function lin.anc.ts() can be applied to data. By setting degree = 0, ancestor regression for i.i.d. linear structural equation models is covered as well <br>
helpers.R: contains function for downstream analysis after getting the p-values from ancestor regression <br>

To create the data for the figures in the paper, run the function stored in anc_rand_simulation.R with its default parameter <br>
The called plotting functions with explaining comments are stored in figures.R. Additional functions to create the ROC are in helpers-figures.R <br>

To recover the data analysis in Section 4, run data_analysis.R in data/. Make sure to store the data in the proper location before. <br>

la-example.R contains a simple example of a simulation with SVAR and prints some summary statistics. Adapt parameters such as the sample size or the noise distributions to analyze various scenarios.
