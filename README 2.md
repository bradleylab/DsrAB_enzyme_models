Sulfur isotope fractionation by dissimilatory sulfite reductase

William D. Leavitt, Alexander S. Bradley, André Santos, Inês A.C. Pereira, David T. Johnston.

These files are the R code files used to process the isotope data from enzyme experiments and solve for fractionation factors (alpha) as well as the 33.lambda exponent relating 34S and 33S fractionation. 

Derivation of the equations is given in the supplementary files to the manuscript. 

Questions regarding these files should be addressed to

Alex Bradley abradley@eps.wustl.edu
or
Wil Leavitt wleavitt@fas.harvard.edu

The included files are:

_Code_for_figures.R_ : this file was used to generate Figure 2

_DsrAB_Model1.R_ : this file runs the calculations to solve for the fractionation between sulfite and sulfonate, and between sulfite and reduced S. It also solves for lambdas. 

Dsr_data_23Aug.txt: a dataframe containing the raw isotope and concentration data

MixingModel.R: this file runs the calculations to solve for the secondary fractionation in the two-fractionation model described in the manuscript supplement (and plotted in Figure 2C)

README.txt: this file