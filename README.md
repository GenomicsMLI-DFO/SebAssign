# SebAssign

Pipeline developed to provide species diagnostic for *Sebastes mentella* and S. fasciatus from four microsatellite loci

This pipeline was developed in 2018 around the `assignPOP` R package (Chen et al. 2018) which perferm genetic assignment using machine learning algorithm.

## How to install SebAssign

1. Install R and Rstudio

2. Clone or download this repository:

`https://github.com/GenomicsMLI-DFO/SebAssign`

3. Install the depending R package : `assignPOP` , `readxl` , `dplyr` , `here`, `magrittr`, `tidyr`, `stringr`

```{r}
install.packages(c("assignPOP", "readxl", "dplyr", "here", "magrittr", "tidyr", "stringr"))
```

## How to use it

1. Open the pipeline as an Rstudio project

2. Prepare your files:

-   A assignment files (individuals from unknown population), with the columns **ID** and **2 columns by locus** (with ".A" and ".B" in the column name). Put this file in the folder **02_Data_to_Assign**.

In all cases, make sure that allele character length is always 3 (if it's not the case some minor adjustment are needed)

3.  Set the pipeline

-   Open the **EasyAssign.R** file (preferably within the project in Rstudio).
-   Set the *assign.excel* variable in the parameters section.

4.  Run the pipeline!

## About the reference samples

Reference samples are in the file Ref_Mentella_Fasciatus.gen. They are composed of 100 fasciatus and 81 mentella gulf identified as pure genotype from RADseq redfish dataset (Nesmestan et al. 2020), and genotyped at 4 microsatellite loci ("SEB25", "SEB31", "SEB33", "SEB9", Roques et al. 1999).

This file can be updated if necessary from an excel file (similar to the one describe above) with the function *Update_ref*.

Theses reference samples were validated HOW

## References

Benestan, L., Rougemont, Q., Senay, C., Normandeau, E., Parent, E., Rideout, R., Bernatchez, L., Lambert, Y., Audet, C., and Parent, G.J. 2020. Population genomics and history of speciation reveal fishery management gaps in two related redfish species (*Sebastes mentella* and *Sebastes fasciatus*). Evol. Appl. https://doi.org/10.1111/eva.13143

Chen, K.-Y., Marschall, E.A., Sovic, M.G., Fries, A.C., Gibbs, H.L., and Ludsin, S.A. 2018. assignPOP: An R package for population assignment using genetic, non-genetic, or integrated data in a machine-learning framework. Methods Ecol. Evol. 9, 439-446.

Roques, S., Pallotta, D., Sévigny, J.-M., and Bernatchez, L. 1999a. Isolation and characterisation of microsatellite markers in the North Atlantic redfish (Teleostei : Scorpaenidae, genus *Sebastes*). Mol. Ecol. 8, 685-702.
