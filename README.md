# SebAssign

Pipeline developed to provide species diagnostic for *Sebastes mentella* and S. fasciatus from four microsatellite loci

This pipeline was developed to XYZ.

## Before running the pipeline

You will need to install a few packages: `assignPOP` , `readxl` , `dplyr` , `here`, `magrittr`, `tidyr`, `stringr`

```{r}

install.packages(c("assignPOP", "readxl", "dplyr", "here", "magrittr", "tidyr", "stringr"))

```

## How to use it

1.  Clone or download this repository:

`https://github.com/GenomicsMLI-DFO/SebAssign`

2.  Prepare your files:

-   A assignment files (individuals from unknown population), with the columns **ID** and **2 columns by locus** (with ".A" and ".B" in the column name). Put this file in the folder **02_Data_to_Assign**.

In all cases, make sure that allele character length is always 3 (if it's not the case some minor adjustment are needed)

3.  Set the pipeline

-   Open the **EasyAssign.R** file (preferably within the project in Rstudio).
-   Set the *assign.excel* variable in the parameters section.

4.  Run the pipeline!

## About the reference samples

Reference samples are in the file Ref_Mentella_Fasciatus.gen. They are composed of 100 fasciatus and 81 mentella gulf identified as pure from the RADseq redfish papers, and genotyped at 4 microsatellites loci ("SEB25", "SEB31", "SEB33", "SEB9").

This file can be updated if necessary from an excel file (similar to the one describe above) with the function *Update_ref*.

Theses reference samples were validated ...
