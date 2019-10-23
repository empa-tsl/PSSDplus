# PSSDplus
Generates probabilistic species sensitivity distributions. 

PSSD+ is a tool developed as part of the following publication:
Systematic consideration of parameter uncertainty and variability in probabilistic species sensitivity distributions
Authors: Henning Wigger, Delphine Kawecki, Bernd Nowack and Véronique Adam
Institute: Empa, Swiss Federal Laboratories for Materials Science and Technology, Technology and Society Laboratory, Lerchenfeldstrasse 5, 9014 St. Gallen, Switzerland
Published in Integrated Environmental Assessment and Management in 2019. 

This work is licensed under: Creative Commons Attribution Non Commercial Share Alike 4.0 International

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3516566.svg)](https://doi.org/10.5281/zenodo.3516566)

The PSSD+ tool generates probabilistic species sensitivity distributions (PSSDs). 
It includes five functions: 
- rmore creates no-observed effect concentration (NOEC) distributions combining NOEC values for species for which 3 or more data points are available. This function is used by the functions do.pSSD, do.pSSD.Ag, do.pSSD.troph and do.pSSD.troph.Ag. 
- do.pSSD builds the PSSD in a generic case, where coefficients of variations are defined for each uncertainty factor.
- do.pSSD.troph builds the PSSD for a weighted dataset, representing a typical structure of freshwater communities.
- do.pSSD.Ag and do.pSSD.Ag.troph are equivalent to do.pSSD and do.pSSD.troph, respectively, but they were built to work in the specific case of Ag, where the distribution associated with the uncertainty factor used to convert acute data to chronic values is not defined by a modal value and a coefficient of variation. Instead, a minimum value of 1, a mode of 10 and a maximum value of 100 were used.

The functions do.pSSD, do.pSSD.troph, do.pSSD.Ag and do.pSSD.Ag.troph take as inputs three csv files:
- A table of ecotoxicity values (data points) extracted from the literature
- A table of uncertainty factors associated with these values to correct for the exposure duration (〖UF〗_t). These are defined following the methodology detailed in Wigger et al.
- A table of uncertainty factors associated with the ecotoxicity values to convert them to NOECs (〖UF〗_dd). These are defined following the methodology detailed in Wigger et al.
In these files, the first row gives the names of the species. One column represents one species and the data points or UFs are given on the rows. Each cell in these three files represents the same data point. 
An additional csv file is needed for do.pSSD.troph and do.pSSD.Ag.troph, where the trophic level of each species is represented by a number. 1 represents primary producers, 2 represents herbivores and 3 represents carnivores. The first row

do.pSSD and do.pSSD.Ag follow three main calculation steps:
1. Data points retrieved from the literature are converted to chronic NOEC values, using uncertainty factors.
2. For each species, a probability distribution is built that represents the NOEC values and their associated uncertainties. When one data point is available for a species, a triangular distribution is built. When two data points are available for a species, a trapezoidal distribution is built. When three or more data points are available for one species, a more complex distribution is built, using the function rmore.
3. All species-specific distributions are combined in a pSSD. This step is repeated SIM times.

do.pSSD.troph and do.pSSD.Ag.troph follow the same calculation steps, to which a weighting step is added. This additional step allows sampling species so that typical proportions of primary producers, herbivores and carnivores are represented on the pSSD. 
