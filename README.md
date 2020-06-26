# polyan: a python package for modelling polysome profiles from ribosome density data

## Purpose

polyan is designed to produce visual output resembling polysome profiles generated with sucrose density gradients. 
These modelled profiles can serve as a convenient quality checking tool for data relating to ribosome densities on mRNAs, and can be used to generate secondary visualisations for assessing movement trends of transcripts between polysome peaks. 

## Data requirements

polyan takes as input pandas dataframes which contain, as minimum requirements, a column specifying transcript names, and a second column specifying relative ribosome densities on each transcript. The most common source of input data for this analysis would likely ribosome footprinting experiemnts, in which case the second column should specify read counts for each transcript (and importantly, shoudl not be RPKM). Where a minimum dataset with only transcript names and ribosome densities is used, transcript abundance data are added from generic datasets provided by polyan. These are currently available for *Saccharomyces cerevisiae* cells and for HEK293 cells. ALternatively, the dataframe can contain a third column specifying transcript abundances (if these are from sequencing data they should again be in the form of read counts rather than RPKM).

## Tutorials

Detailed tutorials are in preparation. 