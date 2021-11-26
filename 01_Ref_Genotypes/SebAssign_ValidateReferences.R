# Info --------------------------------------------------------------------

# Pipeline developed to provide species assignement for Sebastes mentella 
# and S. fasciatus, using 4 microsatelittes loci 
# 
# THIS SPECIFIC PART ALLOW TO TEST THE REFERENCE SAMPLES
# It only need to be run one time
#
# Writen by: Audrey Bourret
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
for(i in 1:length( list.files("./03_Functions") )){
  source(file.path("./03_Functions",  list.files("./03_Functions")[i]))  
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



# Table 1 -----------------------------------------------------------------


res.files <- list.files(file.path(here::here(), "01_Ref_Genotypes", "Ref_Mentella_Fasciatus", "MC_cross-validation_naiveBayes"), pattern = "Out_0.._1", full.names = T)

res <- data.frame()

for(x in res.files){
 
  int <- read.delim(res.files[1], sep="")
 
  res <- bind_rows(res, int)
  
}

head(res)
n.tests <- nrow(res)

res <- res %>% mutate(pred.pop.50 = ifelse(between(pop.1, 0.5, 0.5), "undetermined", pred.pop),
               pred.pop.55 = ifelse(between(pop.1, 0.45, 0.55), "undetermined", pred.pop),
               pred.pop.60 = ifelse(between(pop.1, 0.40, 0.60), "undetermined", pred.pop),
               pred.pop.65 = ifelse(between(pop.1, 0.35, 0.65), "undetermined", pred.pop),
               pred.pop.70 = ifelse(between(pop.1, 0.30, 0.70), "undetermined", pred.pop),
               pred.pop.75 = ifelse(between(pop.1, 0.25, 0.75), "undetermined", pred.pop),
               pred.pop.80 = ifelse(between(pop.1, 0.20, 0.80), "undetermined", pred.pop),
               pred.pop.85 = ifelse(between(pop.1, 0.15, 0.85), "undetermined", pred.pop),
               pred.pop.90 = ifelse(between(pop.1, 0.10, 0.90), "undetermined", pred.pop),
               pred.pop.95 = ifelse(between(pop.1, 0.05, 0.95), "undetermined", pred.pop)) 

res <- res %>% select(-c("pred.pop", "pop.1", "pop.2")) %>% 
  pivot_longer(cols = names(.) %>% str_subset("Ind.ID|origin.pop", negate = T),
               names_to = "threshold",
               values_to = "pred.pop")
table.1 <- res %>% mutate(Test = ifelse(origin.pop == pred.pop, "TRUE", "FALSE"),
               Assigned = ifelse(pred.pop %in% c("pop.1", "pop.2"), "assigned", pred.pop)) %>% 
        group_by(threshold, Assigned, Test) %>% summarise(N = n()) %>% 
        mutate(Group = ifelse(Assigned == "assigned", Test, Assigned),
               threshold = str_remove(threshold, "pred.pop.")) %>% 
        group_by(threshold, Group) %>% summarise(N = sum(N)/n.tests) %>% 
        pivot_wider(names_from = Group, values_from = N, values_fill = 0) 

table.1

# Print in md

library(knitr)

table1.html <- knitr::kable(table.1, digits = 3, format = "html")

cat(table1.html, file = file.path(here::here(), "01_Ref_Genotypes", "Ref_Mentella_Fasciatus", "Table1.html"))
