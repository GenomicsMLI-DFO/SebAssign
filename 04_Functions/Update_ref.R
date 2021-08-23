# Function to update reference files



update_ref <- function(locus = c("SEB25", "SEB31", "SEB33", "SEB9"),
                       pop.name =  c("fasciatus", "mentella_golf"),
                       ref.excel = file.path(here::here(),"01_Ref_Genotypes", "RAD2019_4groupes_sebastes.xlsx"),
                       ref.file = file.path(here::here(), "01_Ref_Genotypes", "Ref_Mentella_Fasciatus.gen")
                       ){
  
  if(file.exists(ref.excel)) stop("There is no Excel file within the 01_Ref_Genotype for the update")
  
  # Step 1 : Excel to Dataframe format
  ref.geno <- base::as.data.frame(readxl::read_excel(ref.excel,  col_types = rep("text", 2 + length(locus) * 2)))
  
  # Keep only the selected population
  ref.geno <- ref.geno %>% dplyr::filter(POP %in% pop.name)
  
  ref.df <- merge.MSAT.alleles(data=ref.geno, locus=locus, na="NA")
  ref.df
  
  # Step 2: Save as genpop
  
  write.genpop(fn = ref.file, 
               data = ref.df, 
               pop = "POP", 
               ind = "ID",
               locus = locus)
    
  
}
