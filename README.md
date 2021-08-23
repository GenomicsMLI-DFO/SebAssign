# SebAssign
Pipeline developed to provide species diagnostic for Sebastes mentella and S. fasciatus

## How to use it

1. Clone or download this repository:

`https://github.com/biodray/EasyAssign.git`

2. Prepare your files:

- A reference file in Excel, with the columns **POP**, **ID**, then **2 columns by locus** (with ".A" and ".B" in the column name). Put this file in the folder **01_Ref_Genotypes**.
- A assignment files (individuals from unknown population), with the columns **ID** and **2 columns by locus** (with ".A" and ".B" in the column name). Put this file in the folder **02_Data_to_Assign**.

In all cases, make sure that allele character length is always 3 (if it's not the case some minor adjustment are needed)

3. Prepare Rstudio

- Open the **EasyAssign.Rproj** file in Rstudio (to work in right working directory).
- Open the **EasyAssign.R** file in Rstudio.
- Set the *locus*, *ref.excel* and *assign.excel* variable in the parameters section.

4. Run the pipeline!


## About the reference samples

Reference samples are in the file Ref_Mentella_Fasciatus.gen. They are composed of 100 fasciatus and 81 mentella gulf identified as pure from the RADseq redfish papers, and genotyped at 4 microsatellites loci ("SEB25", "SEB31", "SEB33", "SEB9"). 

This file can be updated if necessary from an excel file (similar to the one describe above) with the function *Update_ref*. 

Theses reference samples were validated ...