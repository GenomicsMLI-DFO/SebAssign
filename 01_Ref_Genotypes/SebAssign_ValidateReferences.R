# Info --------------------------------------------------------------------

# Pipeline developed to provide species diagnostic for Sebastes mentella 
# and S. fasciatus, using 4 microsatelittes loci 
# 
# THIS SPECIFIC PART ALLOW TO TEST THE REFERENCE SAMPLES
# It only need to be run one time
#
# Write by: Audrey Bourret
#

# See https://alexkychen.github.io/assignPOP/index.html for more information
# on the assignPOP package

# Library -----------------------------------------------------------------

# install.packages(c("assignPOP", "readxl", "dplyr", "here", "magrittr", "tidyr", "stringr"))

# POPassign 
library(assignPOP)

# Tools to play around files
library(readxl)
library(dplyr)
library(tidyr)
library(here)
library(magrittr)
library(stringr)

#if(!require(klaR)){ install.packages("klaR") }
#library(klaR)

# Load internal functions
for(i in 1:length( list.files("./04_Functions") )){
  source(file.path("./04_Functions",  list.files("./04_Functions")[i]))  
}


# Parameters --------------------------------------------------------------

# Set these parameters

locus <- c("SEB25", "SEB31", "SEB33", "SEB9")
ref.gen    <- assignPOP::read.Genepop(file.path(here::here(), "01_Ref_Genotypes", "Ref_Mentella_Fasciatus.gen"), pop.names= c("fasciatus", "mentella_golf"), haploid = FALSE)

# Evaluate baseline -------------------------------------------------------

# Compute cross-validation statistics

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

