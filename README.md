# SebAssign

Pipeline developed to provide species assignment for *Sebastes mentella* and *S. fasciatus* from four microsatellite loci


__Main author:__  Audrey Bourret    
__Affiliation:__  Fisheries and Oceans Canada (DFO)   
__Group:__        Laboratory of genomics   
__Location:__     Maurice Lamontagne Institute  
__Affiliated publication:__  [Senay, C., Bermingham, T., Parent, G.J., Benoît, H. P., Parent, E., Bourret, A. 2022. Identifying two Redfish species, Sebastes mentella and S. fasciatus, in fishery and surveycatches using anal fin ray count in Units 1 and 2. Can. Tech. Rep. Fish. Aquat. Sci. 3445: viii +46 p.](https://waves-vagues.dfo-mpo.gc.ca/Library/41043364.pdf) 
__Contact:__      audrey.bourret@dfo-mpo.gc.ca

- [Objective](#objective)
- [Status](#status)
- [Versions](#versions)
- [Requirements](#requirements)
- [How to use SebAssign](#how-to-use-sebassign)
- [Reference samples](#about-the-reference-samples)
- [References](#references)

## Objective

This pipeline was developed in 2019, and make use of the `assignPOP` R package (Chen et al. 2018) alowing users to perform genetic assignment using diverse machine learning algorithms. The ultimate goal of this pipeline was to ease and standardise the genetic identification for stock management of these species. 

## Status
Completed

## Versions

Version 0.1.0 can be uploaded [here](https://github.com/GenomicsMLI-DFO/SebAssign/releases/tag/v0.1.0)

Development version can be uploaded directly using the **Code** green buttom above!

## Requirements

1. Install R and Rstudio

2. Clone or download this repository: https://github.com/GenomicsMLI-DFO/SebAssign

3. Install the depending R package : `assignPOP` , `readxl` , `dplyr` , `here`, `magrittr`, `tidyr`, `stringr`, `ggplot2`

This can be done all at once with this command line in R :

```{r}
install.packages(c("assignPOP", "readxl", "dplyr", "here", "magrittr", "tidyr", "stringr", "ggplot2"))
```

## How to use SebAssign

1. Open the pipeline as an Rstudio project

2. Prepare your files:

-   A assignment files (individuals from unknown population), with the columns **ID** and **2 columns by locus** (with ".A" and ".B" in the column name). Put this file in the folder **02_Data_to_Assign**. The pipeline was made to read excel files directly (using the `readxl` package). 

In all cases, make sure that allele character length is always 3 (e.g. do not code 99 as allele, but 099)

3.  Set the pipeline

-   Open the **EasyAssign.R** file (preferably within the project in Rstudio).
-   Set the *assign.excel* variable in the parameters section.

4.  Run the pipeline! Results will be save in the **02_Data_to_Assign** folder, in a sub directory named after the excel file. 

A test dataset (*RAD2019_sebastes_Test.xlsx*) is included with the pipeline, as well as the its assignment results.

![figTest](/02_Data_to_Assign/RAD2019_sebastes_Test/Results_assignments.png)


## About the reference samples

Reference samples are in the file **Ref_Mentella_Fasciatus.gen**. They are composed of 100 fasciatus and 81 mentella gulf identified as pure genotype from RADseq redfish dataset (Q-value > 0.95 in Benestan et al. 2020), and genotyped at 4 microsatellite loci (*SEB9*, *SEB25*, *SEB31*, *SEB33*; Roques et al. 1999).

This file can be updated if necessary from an excel file (similar to the one describe above) with the function *Update_ref*.

Theses reference samples were tested under five different machine-learning algorithm (LDA, Naive Bayes, Random Forest, SVM and Decision tree) using MCMC cross-validation tests (training set : 0.5, 0.7 or 0.9 of all reference samples, 100 iterations). Both the Naive Bayes and the SVM model showed best results, and **Naive Bayes** was chosen to perform genetic assignment of samples from unknown species.  

![fig1](/01_Ref_Genotypes/Ref_validation_MCMC.png)

From MCMC cross-validation results using Naive Bayes model, we tested the impact of increasing the membership probability threshold (from 50% to 95%) on the assignment rates.  

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> threshold </th>
   <th style="text-align:right;"> FALSE </th>
   <th style="text-align:right;"> TRUE </th>
   <th style="text-align:right;"> undetermined </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 50 </td>
   <td style="text-align:right;"> 0.077 </td>
   <td style="text-align:right;"> 0.923 </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 55 </td>
   <td style="text-align:right;"> 0.066 </td>
   <td style="text-align:right;"> 0.901 </td>
   <td style="text-align:right;"> 0.033 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 60 </td>
   <td style="text-align:right;"> 0.055 </td>
   <td style="text-align:right;"> 0.879 </td>
   <td style="text-align:right;"> 0.066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 65 </td>
   <td style="text-align:right;"> 0.055 </td>
   <td style="text-align:right;"> 0.868 </td>
   <td style="text-align:right;"> 0.077 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 70 </td>
   <td style="text-align:right;"> 0.055 </td>
   <td style="text-align:right;"> 0.868 </td>
   <td style="text-align:right;"> 0.077 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 75 </td>
   <td style="text-align:right;"> 0.055 </td>
   <td style="text-align:right;"> 0.868 </td>
   <td style="text-align:right;"> 0.077 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 80 </td>
   <td style="text-align:right;"> 0.044 </td>
   <td style="text-align:right;"> 0.857 </td>
   <td style="text-align:right;"> 0.099 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 85 </td>
   <td style="text-align:right;"> 0.033 </td>
   <td style="text-align:right;"> 0.824 </td>
   <td style="text-align:right;"> 0.143 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 90 </td>
   <td style="text-align:right;"> 0.022 </td>
   <td style="text-align:right;"> 0.813 </td>
   <td style="text-align:right;"> 0.165 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 95 </td>
   <td style="text-align:right;"> 0.011 </td>
   <td style="text-align:right;"> 0.747 </td>
   <td style="text-align:right;"> 0.242 </td>
  </tr>
</tbody>
</table>

The code to perform these analysis can be find here: **/01_Ref_Genotypes/SebAssign_ValidateReferences.R**

## References

Benestan, L., Rougemont, Q., Senay, C., Normandeau, E., Parent, E., Rideout, R., Bernatchez, L., Lambert, Y., Audet, C., and Parent, G.J. 2021. Population genomics and history of speciation reveal fishery management gaps in two related redfish species (*Sebastes mentella* and *Sebastes fasciatus*). Evol. Appl. https://doi.org/10.1111/eva.13143

Chen, K.-Y., Marschall, E.A., Sovic, M.G., Fries, A.C., Gibbs, H.L., and Ludsin, S.A. 2018. assignPOP: An R package for population assignment using genetic, non-genetic, or integrated data in a machine-learning framework. Methods Ecol. Evol. 9, 439-446.

Roques, S., Pallotta, D., Sévigny, J.-M., and Bernatchez, L. 1999a. Isolation and characterisation of microsatellite markers in the North Atlantic redfish (Teleostei : Scorpaenidae, genus *Sebastes*). Mol. Ecol. 8, 685-702.
