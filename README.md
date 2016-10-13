# DrosophilaWGCNA
R scripts used to analyse data from Innocenti and Morrow 2010 using WGCNA, and more!

This repo relates to an in-prep project analysing data from [Innocenti and Morrow 2010 PLoS Biology](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000335). The analysis makes heavy use of the WGCNA package ([see tutorials here](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/)).

The file "allData.R" is an R save file containing an R list with multiple files needed to reproduce the analyses and figures. This file was created using the "Innocenti networks - make data.R" script, which in turn relies on functions found in the "functions.R" script and the "Normalise_expression_data.R" script. The "Innocenti networks - data analysis.R" script can be used to generate all the statistical analyses and figures in the forthcoming paper, and relies upon the "functions.R" and "allData.R" files.

The R scripts are heavily annotated, so I've archived them here in case they are useful to anyone learning to use WGCNA.
