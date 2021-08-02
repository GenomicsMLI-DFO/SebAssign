# Info --------------------------------------------------------------------

# Perform assignment with assignPOP
# 
# Write by: Audrey Bourret
# 2019-10-15
#

# See https://alexkychen.github.io/assignPOP/index.html

# Library -----------------------------------------------------------------

# POPassign - ** Attention R >= 3.6.0 **
if(!require(assignPOP)){ install.packages("assignPOP") }
library(assignPOP)

if(!require(klaR)){ install.packages("klaR") }
library(klaR)

# readxl convert read gen data
if(!require(readxl)){ install.packages("readxl") }
library(readxl)

#tidyverse to play around xl filesconvert read gen data
if(!require(tidyverse)){ install.packages("tidyverse") }
library(tidyverse)

# Internal functions
for(i in 1:length( list.files("./04_Functions") )){
  source(file.path("./04_Functions",  list.files("./04_Functions")[i]))  
}


# Parameters --------------------------------------------------------------

# Set these parameters

locus <- c("SEB25", "SEB31", "SEB33", "SEB9")

#ref.excel <- "STD_Assignation_men-Fasc_3130_EP_14052019.xlsx" # Vieilles séquences références
ref.excel <- "RAD2019_4groupes_sebastes.xlsx" # Séquences références validées avec les données RAD

assign.excel <- "TEL15_Coh16_BIN-MUX1-ASSIG_final_2019.xlsx"

# Nothing to change here
ref.file  <- ref.excel %>% str_replace(".xlsx", ".gen") %>% str_replace(".xls", ".gen")
ref.dir   <- file.path("01_Ref_Genotypes" ,ref.excel %>% str_remove(".xlsx") %>% str_remove(".xls"))

assign.file  <- assign.excel %>% str_replace(".xlsx", ".gen") %>% str_replace(".xls", ".gen")
assign.dir   <- file.path("03_Results" ,assign.excel %>% str_remove(".xlsx") %>% str_remove(".xls"))

# Prepare Genetic Data ----------------------------------------------------

# Step 1 : Excel to Dataframe format
ref.geno <- as.data.frame(read_excel(file.path("01_Ref_Genotypes",ref.excel),  col_types = rep("text", 2 + length(locus) * 2)))

# Keep only 2 populations
ref.geno <- ref.geno %>% filter(POP %in% c("fasciatus", "mentella_golf"))

ref.df <- merge.MSAT.alleles(data=ref.geno, locus=locus, na="NA")
ref.df


# Step 2 : Merge alleles
assign.geno <- as.data.frame(read_excel(file.path("02_Data_to_Assign",assign.excel),  col_types = rep("text", 2 + length(locus) * 2)))

assign.geno <- assign.geno %>% drop_na()

assign.df <- merge.MSAT.alleles(data=assign.geno, locus=locus, na="NA")
assign.df

# Step 3: Save as genpop
write.genpop(fn = file.path("01_Ref_Genotypes", ref.file), 
             data = ref.df, 
             pop = "POP", 
             ind = "ID",
             locus = locus)

write.genpop(fn = file.path("02_Data_to_Assign",assign.file), 
             data = assign.df, 
             ind = "ID",
             locus = locus)

# Step 4: Upload in the rigth format

ref.gen    <- read.Genepop( file.path("01_Ref_Genotypes",ref.file), pop.names=unique(ref.df$SP), haploid = FALSE)
assign.gen <- read.Genepop( file.path("02_Data_to_Assign",assign.file), pop.names="POP1", haploid = FALSE)

# Evaluate baseline -------------------------------------------------------

# remove low variance loci

#ref.gen.rd <- reduce.allele(ref.gen, p = 0.95)

# NOTE: it is not necessary to run this part each time ...

# Compute cross-validation statistics

# Make a directory
dir.create(ref.dir)


# Population assignment test using Monte-Carlo cross-validation

#"lda", "svm", "naiveBayes", "tree", and "randomForest"

assign.MC(ref.gen, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.75, 1),
           loci.sample="fst", iterations=100, model="svm", dir=paste0(file.path(ref.dir,"MC_cross-validation_svm"),"/"))

assign.MC(ref.gen, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.75, 1),
          loci.sample="fst", iterations=100, model="naiveBayes", dir=paste0(file.path(ref.dir,"MC_cross-validation_naiveBayes"),"/"))

assign.MC(ref.gen, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.75, 1),
          loci.sample="fst", iterations=100, model="tree", dir=paste0(file.path(ref.dir,"MC_cross-validation_tree"),"/"))

assign.MC(ref.gen, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.75, 1),
          loci.sample="fst", iterations=100, model="randomForest", dir=paste0(file.path(ref.dir,"MC_cross-validation_randomForest"),"/"))


accuMC.svm <- accuracy.MC(dir = paste0(file.path(ref.dir,"MC_cross-validation_svm"),"/"))
accuMC.naiveBayes <- accuracy.MC(dir = paste0(file.path(ref.dir,"MC_cross-validation_naiveBayes"),"/"))
accuMC.tree <- accuracy.MC(dir = paste0(file.path(ref.dir,"MC_cross-validation_tree"),"/"))
accuMC.randomForest <- accuracy.MC(dir = paste0(file.path(ref.dir,"MC_cross-validation_randomForest"),"/"))

accuracy.plot(accuMC.svm, pop = c("all", "pop.1", "pop.2")) +
   geom_hline(yintercept = 0.9, lty = "dashed", col = "red") +
   geom_hline(yintercept = 0.95, lty = "dashed", col = "blue") +
   ylim(0.5, 1) +
   ggtitle("SVM")

accuracy.plot(accuMC.naiveBayes, pop = c("all", "pop.1", "pop.2")) +
   geom_hline(yintercept = 0.9, lty = "dashed", col = "red") +
   geom_hline(yintercept = 0.95, lty = "dashed", col = "blue") +
   ylim(0.5, 1) +
   ggtitle("naiveBayes")
   
accuracy.plot(accuMC.tree, pop = c("all", "pop.1", "pop.2")) +
   geom_hline(yintercept = 0.9, lty = "dashed", col = "red") +
   geom_hline(yintercept = 0.95, lty = "dashed", col = "blue") +
   ylim(0.5, 1) +
   ggtitle("tree")

accuracy.plot(accuMC.randomForest, pop = c("all", "pop.1", "pop.2")) +
   geom_hline(yintercept = 0.9, lty = "dashed", col = "red") +
   geom_hline(yintercept = 0.95, lty = "dashed", col = "blue") +
   ylim(0.5, 1) +
   ggtitle("randomForest")

# Population assignment test using K-fold cross-validation

assign.kfold(ref.gen, k.fold=c(3, 4, 5), train.loci=c(1), 
             loci.sample="random", model="naiveBayes", dir=paste0(file.path(ref.dir,"kfold_cross-validation"),"/"))

accuKF <- accuracy.kfold(dir = paste0(file.path(ref.dir,"kfold_cross-validation"),"/"))

accuracy.plot(accuKF, pop = c("all", "pop.1", "pop.2")) +
   geom_hline(yintercept = 0.9, lty = "dashed", col = "red") +
   geom_hline(yintercept = 0.96, lty = "dashed", col = "blue") +
   ylim(0.6, 1) +
   ggtitle("naiveBayes")


# Predict sources of unknown individuals ----------------------------------

# Make a directory
dir.create(assign.dir)

# 1.Perform assignment test using genetic data and naive Bayes
assign.X( x1=ref.gen, x2=assign.gen, dir=paste0(file.path(assign.dir,"naiveBayes"),"/"), model="naiveBayes", mplot = T)

# 2.Perform assignment test using decision tree
assign.X( x1=ref.gen, x2=assign.gen, dir=paste0(file.path(assign.dir,"tree"),"/"),  model="tree", mplot = T)

# 3.Perform assignment test using random forest
assign.X( x1=ref.gen, x2=assign.gen, dir=paste0(file.path(assign.dir,"randomForest"),"/"),  model="randomForest", mplot = T)

# 4.Perform assignment test using SVM
assign.X( x1=ref.gen, x2=assign.gen, dir=paste0(file.path(assign.dir,"svm"),"/"),  model="svm", mplot = T)


membership.plot(dir = paste0(file.path(assign.dir,"tree"),"/"))


# Create final file -------------------------------------------------------

res1 <- read.table(file.path(assign.dir,"naiveBayes", "AssignmentResult.txt"), header = T)
res1$Model <- "naiveBayes"

res2 <- read.table(file.path(assign.dir,"tree", "AssignmentResult.txt"), header = T)
res2$Model <- "tree"

res3 <- read.table(file.path(assign.dir,"randomForest", "AssignmentResult.txt"), header = T)
res3$Model <- "randomForest"

res4 <- read.table(file.path(assign.dir,"svm", "AssignmentResult.txt"), header = T)
res4$Model <- "svm"

res <- rbind(res1, res2, res3, res4)


res %>% ggplot(aes(x=pop.1)) + 
               geom_histogram() + 
               geom_vline(xintercept = c(0.05, 0.95), col = "blue") +
               geom_vline(xintercept = c(0.10, 0.90), col = "red") +
   facet_wrap(~Model) +            
   theme_bw()


names(res) <- c("ID", "POP.assign", "fasciatus", "mentella_golf", "model")


res %>% filter(model == "naiveBayes") %>% 
        gather("fasciatus", "mentella_golf", key = "Species", value = "Prob") %>% 
        left_join(assign.df %>% select(ID, POP)) %>% 
        filter(Prob >= 0.97) %>% 
        group_by(POP, model, Species) %>% 
        summarise(N = n()) %>% View()

library(corrplot)

res %>% select(-POP.assign) %>% 
   gather("fasciatus", "mentella_golf", key = "Species", value = "Prob") %>% 
   spread(model, value = Prob) %>%  
    left_join( assign.df %>% select(ID, POP)) %>% 
   
       filter(Species  == "fasciatus",
              POP == "Maria",
              naiveBayes >=0.97) %>% # View()
   select(naiveBayes, randomForest, svm, tree) %>% 
   plot()
   


names(res1) <- c("ID", "POP.assign", "fasciatus", "mentella_golf", "model")

res.final <- res1 %>% mutate(POP.assign = ifelse(POP.assign == "pop.1", "fasciatus", "mentella_golf"),
                                                 POP = ifelse(str_detect(ID, "Maria"), "Maria", "Repro"))
res.final

res.Maria <- res.final %>% filter(POP == "Maria",
                                  mentella_golf >= 0.95) %>% 
                           arrange(desc(mentella_golf)) %>% 
                           transmute(ID, SP=POP.assign, Membership.prob = mentella_golf)
res.Maria    

write.csv2(res.Maria, "Resultats_Mentella_Maria_2019-06-06.csv", row.names = F)

res.Maria2 <- res.final %>% filter(POP == "Maria") %>% 
                            mutate(SP.95 = ifelse(POP.assign == "fasciatus" & fasciatus >= 0.95, "fasciatus",
                                                  ifelse(POP.assign == "mentella_golf" & mentella_golf >= 0.95, "mentella_golf",
                                                         "undetermined")) ) %>% 
                            transmute(ID, SP.95, Prob.fasciatus = fasciatus, Prob.mentella_golf = mentella_golf)

all.ID.Maria <- as.data.frame(read_excel(file.path("02_Data_to_Assign",assign.excel),  col_types = rep("text", 2 + length(locus) * 2))) %>% 
                            select(ID) %>% 
                            filter(str_detect(ID, "Maria"))

res.Maria2 <- res.Maria2 %>% full_join(all.ID.Maria) %>% mutate(SP.95 = ifelse(is.na(SP.95),"undetermined", SP.95))    
res.Maria2

write.csv2(res.Maria2, "Resultats_Mentella_Maria_Complet_2019-06-11.csv", row.names = F)




res.Caro <- res.final %>% filter(POP != "Maria") %>% 
                          mutate(SP.95 = ifelse(POP.assign == "fasciatus" & fasciatus >= 0.95, "fasciatus",
                                         ifelse(POP.assign == "mentella_golf" & mentella_golf >= 0.95, "mentella_golf",
                                          "undetermined")) ) %>% 
                          transmute(ID, SP.95, Prob.fasciatus = fasciatus, Prob.mentella_golf = mentella_golf)

res.Caro                                  
write.csv2(res.Caro, "Resultats_Sebastes_Caro_2019-06-06.csv", row.names = F)

