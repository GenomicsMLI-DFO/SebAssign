# Functions useful to convert xls files to others formats

merge.MSAT.alleles <- function(data, locus, allele = c(".A",".B"), name = "", na){
   for (x in 1:length(locus)){
      LOC   <- locus[x]
      COL   <- c(paste0(LOC, allele[1]), paste0(LOC, allele[2]))      
      NAME  <- paste0(locus[x],name)
      LOCUS <- paste0(data[,COL[1]],data[,COL[2]])
      data[,NAME] <- LOCUS 
      data[which(data[,NAME] == paste0(na,na)),NAME] <- "000000"
      data[which(data[,NAME] == 0),NAME] <- "000000"
      
      # Drop old variable names
      data <- data[ , !(names(data) %in% paste0(NAME, allele))]      
      }

   return(data)
   }


write.genpop <- function(fn, data, pop = NULL, ind, locus){

   # header
   cat("Genpop file format",
       locus,
       append = FALSE,
       file = fn,
       sep = "\n")
   
   
   if(!is.null(pop)){
   # loop for each pop
   
   for(x in unique(data[,pop])){
       cat("POP",
           append = TRUE,
           file = fn,
           sep = "\n")
      
       data.int <- data[which(data[, pop] == x),]
      
       # loop for each line
       
       for(y in 1:nrow(data.int)){
           cat(paste(paste0(substr(x,1,3), "_", data.int[y, ind]),
                     ",",
                     paste(data.int[y, locus], collapse = " "),
                     sep = " "),
               append = TRUE,
               file = fn,
               sep = "\n")
            
       }
       
   }
      
   }

   if(is.null(pop)){
   # loop for each ind
       cat("POP",
           append = TRUE,
           file = fn,
           sep = "\n")
      
       data.int <- data
      
       # loop for each line
       
       for(y in 1:nrow(data.int)){
           cat(paste(data.int[y, ind],
                     ",",
                     paste(data.int[y, locus], collapse = " "),
                     sep = " "),
               append = TRUE,
               file = fn,
               sep = "\n")
            
       }
       
   }
   
   
   
}
