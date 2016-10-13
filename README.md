# DrosophilaWGCNA
R scripts used to analyse data from Innocenti and Morrow 2010 using WGCNA, and more!

This repo relates to an in-prep project analysing data from [Innocenti and Morrow 2010 PLoS Biology](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000335). The analysis makes heavy use of the WGCNA package ([see tutorials here](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/)).

The file "allData.R" is an R save file containing an R list with multiple dataframes needed to reproduce the analyses and figures. This file was created using the "Innocenti networks - make data.R" script, which in turn relies on functions found in the "functions.R" script and the "Normalise_expression_data.R" script, as well as the "sampleIDs.csv" data file and the raw microarray data from Innocenti and Morrow (download it from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17013)). The "Innocenti networks - data analysis.R" script can be used to generate all the statistical analyses and figures in the forthcoming paper, and relies upon the "functions.R" and "allData.R" files. The four .csv files with names beginning "transcripts associated..." or "sexually antagonistic transcripts" are suppplementary tables from [Innocenti and Morrow](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000335).

If you'd like to look in detail at the module data, or want to get a list of sexually antagonisitic D. melanogaster genes, I recommend you open up the "allData.R" file using the opening lines of code in "Innocenti networks - data analysis.R".

The R scripts are heavily annotated, and could be useful if you'd like to apply WGCNA to your own datasets.
