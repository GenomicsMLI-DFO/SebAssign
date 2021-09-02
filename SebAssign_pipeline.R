# Info --------------------------------------------------------------------

# Pipeline developed to provide species assignment for Sebastes mentella 
# and S. fasciatus, using 4 microsatelittes loci 
# 
# Write by: Audrey Bourret
#

# See https://alexkychen.github.io/assignPOP/index.html for more information
# on the assignPOP package

# Library -----------------------------------------------------------------

# install.packages(c("assignPOP", "readxl", "dplyr", "here", "magrittr", "tidyr", "stringr", "ggplot2"))

# POPassign 
library(assignPOP)

# Tools to play around files
library(readxl)
library(dplyr)
library(tidyr)
library(here)
library(magrittr)
library(stringr)
library(ggplot2)

# Load internal functions
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
}


# Parameters --------------------------------------------------------------

# Set these parameters

locus <- c("SEB25", "SEB31", "SEB33", "SEB9")
threshold <- 95
   
list.files(file.path(here::here(), "02_Data_to_Assign"))
assign.excel <- "RAD2019_sebastes_Test.xlsx"


# Check if these parameters are OK
# Name of the folder where results will be saved
assign.dir   <- file.path(here::here(), "02_Data_to_Assign", assign.excel %>% stringr::str_remove(".xlsx|.xls"))
assign.dir

assign.file  <- assign.excel %>% stringr::str_replace(".xlsx|.xls", ".gen")
assign.file

# Prepare Genetic Data ----------------------------------------------------

# Step 1 : Excel to Dataframe format

assign.geno <- as.data.frame(readxl::read_excel(file.path(here::here(), "02_Data_to_Assign", assign.excel),  col_types = rep("text", 2 + length(locus) * 2)))
assign.geno <- assign.geno %>% tidyr::drop_na()

# Step 2 : Merge alleles
assign.df <- merge.MSAT.alleles(data=assign.geno, locus=locus, na="NA")
assign.df

# Step 3: Save as genpop

write.genpop(fn = file.path(here::here(), "02_Data_to_Assign", assign.file), 
             data = assign.df, 
             ind = "ID",
             locus = locus)

# Step 4: Upload in the right format

ref.gen    <- assignPOP::read.Genepop(file.path(here::here(), "01_Ref_Genotypes", "Ref_Mentella_Fasciatus.gen"), pop.names= c("fasciatus", "mentella_golf"), haploid = FALSE)

assign.gen <- assignPOP::read.Genepop(file.path(here::here(), "02_Data_to_Assign", assign.file), pop.names="POP1", haploid = FALSE)


# Predict sources of unknown individuals ----------------------------------

# Make a directory
if(file.exists(assign.dir)){
   cat("\nThe folder:", assign.dir, "already exist, nothing was done (be careful to not overwrite previous data)\n")
} else {
  cat("\nThe folder:", assign.dir, "was created")
   dir.create(assign.dir) 
}

# 1.Perform assignment test using genetic data and naive Bayes
assign.X( x1=ref.gen, x2=assign.gen, dir=paste0(file.path(assign.dir,"naiveBayes"),"/"), model="naiveBayes", mplot = T)


# Create final file -------------------------------------------------------

res <- read.table(file.path(assign.dir,"naiveBayes", "AssignmentResult.txt"), header = T)

res <- res %>%  mutate(SP.95 = ifelse(pred.pop == "fasciatus" & fasciatus >= 0.95, "fasciatus",
                               ifelse(pred.pop == "mentella_golf" & mentella_golf >= 0.95, "mentella_golf",
                                          "undetermined")) ) %>%
                left_join(assign.df %>% select(ID, POP), by = c("Ind.ID" = "ID") ) %>% 
                transmute(Ind.ID, POP, SP.95, Prob.fasciatus = fasciatus, Prob.mentella_golf = mentella_golf)

res

write.csv2(res, file.path(assign.dir, paste0("Results_assignments_", Sys.Date(), ".csv")), row.names = F)


# Stats -------------------------------------------------------------------

fig.res <- res %>% group_by(POP, SP.95) %>% 
      summarise(N = n()) %>% 
      ggplot(aes(x = POP, y = N, fill = SP.95)) +
      geom_bar(stat= "identity") + 
      theme_bw()

ggsave(plot = fig.res, filename = file.path(assign.dir, "Results_assignments.png"),
       width = 7,
       height = 4,
       units = c("in"))
