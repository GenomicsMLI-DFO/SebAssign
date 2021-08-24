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
library(dplyr)
library(tidyr)
library(here)
library(magrittr)
library(stringr)
library(ggplot2)

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
ref.dir <- file.path(here::here(), "01_Ref_Genotypes", "Ref_Mentella_Fasciatus")
ref.dir

# Evaluate baseline -------------------------------------------------------

# Compute cross-validation statistics

# Population assignment test using Monte-Carlo cross-validation

# List of available machine-learning algorithm: "lda", "svm", "naiveBayes", "tree", and "randomForest"

assign.MC(ref.gen, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.75, 1),
          loci.sample="fst", iterations=100, model="lda", dir=paste0(file.path(ref.dir,"MC_cross-validation_lda"),"/"))

assign.MC(ref.gen, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.75, 1),
          loci.sample="fst", iterations=100, model="svm", dir=paste0(file.path(ref.dir,"MC_cross-validation_svm"),"/"))

assign.MC(ref.gen, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.75, 1),
          loci.sample="fst", iterations=100, model="naiveBayes", dir=paste0(file.path(ref.dir,"MC_cross-validation_naiveBayes"),"/"))

assign.MC(ref.gen, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.75, 1),
          loci.sample="fst", iterations=100, model="tree", dir=paste0(file.path(ref.dir,"MC_cross-validation_tree"),"/"))

assign.MC(ref.gen, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.75, 1),
          loci.sample="fst", iterations=100, model="randomForest", dir=paste0(file.path(ref.dir,"MC_cross-validation_randomForest"),"/"))


# Create a big ref

accuMC.all <- data.frame()

# Loop over the models
for(m in c("lda", "svm", "naiveBayes", "tree", "randomForest")){
    
    print(m)
    
    accuMC <- accuracy.MC(dir = paste0(ref.dir,"/MC_cross-validation_",m,"/"))
    accuMC <- accuMC %>% dplyr::mutate(model = m)
    
    accuMC.all <- dplyr::bind_rows(accuMC.all, accuMC)
    
  }
  

accuMC.all <- accuMC.all %>% tidyr::pivot_longer(names(accuMC.all) %>% stringr::str_subset("assign.rate"),
                                          names_to = "group", values_to = "assign.rate") %>% 
  dplyr::mutate(group = group %>% stringr::str_remove("assign.rate."),
         assign.rate = as.numeric(as.character(assign.rate)),
          model = as.factor(model)
  ) %>% 
  dplyr::filter(!is.na(assign.rate),
         group != "all"#,
         #str_detect(ref, "2REF")
  ) 


# Change levels of model
levels(accuMC.all$model)
levels(accuMC.all$model) <- c("LDA", "Naive Bayes", "Random Forest", "SVM", "Decision tree")         

# Compute stats

accuMC.all %>% filter(train.loci == "1",
                      train.inds == 0.7) %>% 
  group_by(train.loci, model, group, train.inds) %>% summarise(mean = mean(assign.rate),
                                                                sd = sd(assign.rate))

# Figure

fig1 <- accuMC.all %>% filter(train.loci == "1") %>%  
  ggplot(aes(x = group, y = assign.rate, fill = factor(train.inds))) +
  geom_hline(yintercept = c(0.5), col = "gray75") +

  geom_boxplot() +
  geom_hline(yintercept = c(0.9), lty = "dashed", col = "red") +
  facet_grid(. ~ model, space = "free_x", scale = "free_x") +
  labs(x = "Reference population", y = "Assigment accuracy") +
  scale_x_discrete(labels = c("fasciatus", "mentella")) +
  scale_fill_manual(values = c("gray100", "gray75", "gray40"))+
  guides(fill=guide_legend(title="Training set size")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top")
fig1

ggsave(filename = file.path(here::here(), "01_Ref_Genotypes", "Ref_validation_MCMC.png"), plot = fig1, 
       width = 7, height = 3 , units = "in",
       dpi = 300)


# Population assignment test using K-fold cross-validation

assign.kfold(ref.gen, k.fold=c(3, 4, 5), train.loci=c(1), 
             loci.sample="random", model="naiveBayes", dir=paste0(file.path(ref.dir,"kfold_cross-validation"),"/"))

accuKF <- accuracy.kfold(dir = paste0(file.path(ref.dir,"kfold_cross-validation"),"/"))

accuracy.plot(accuKF, pop = c("all", "pop.1", "pop.2")) +
  geom_hline(yintercept = 0.9, lty = "dashed", col = "red") +
  geom_hline(yintercept = 0.96, lty = "dashed", col = "blue") +
  ylim(0.6, 1) +
  ggtitle("naiveBayes")

