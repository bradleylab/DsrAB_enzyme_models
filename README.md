**Sulfur isotope fractionation by dissimilatory sulfite reductase**

William D. Leavitt, Alexander S. Bradley, André Santos, Inês A.C. Pereira, David T. Johnston.

These files contain the R code used to process the isotope data from enzyme experiments and solve for fractionation factors (alpha) as well as the <sup>33</sup>lambda exponent relating <sup>34</sup>S and <sup>33</sup>S fractionation. 

Derivation of the equations is given in the supplementary files to the manuscript. 

Questions regarding these files should be addressed to

Alex Bradley abradley@eps.wustl.edu
or
Wil Leavitt wleavitt@eps.wustl.edu

The included files are:

***Code\_for\_figures.R***: this file was used to generate Figure 2

***DsrAB\_Model1.R***: this file runs the calculations to solve for the fractionation between sulfite and sulfonate, and between sulfite and reduced S. It also solves for <sup>33</sup>lambdas. 

***Dsr\_data\_23Aug.txt***: a dataframe containing the raw isotope and concentration data

***MixingModel.R***: this file runs the calculations to solve for the secondary fractionation in the two-fractionation model described in the manuscript supplement (and plotted in Figure 2C)

***README.md***: this file
