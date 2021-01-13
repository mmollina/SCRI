require(ggplot2)
require(dplyr)

# Using 'system'
setwd("~/repos/SCRI/supermassa_vcf2sm/ploidy_estimation/")
res<-system("python2.7 ../supermassa/src/SuperMASSA.py --print_genotypes --inference f1 --ploidy_range 6:12 --file sugarcane_data/SugSNP225_progeny --f1_parent_data sugarcane_data/SugSNP225_parents", 
            intern = TRUE)
res

# Windows
# setwd("C:/Users/mmoll/OneDrive/Desktop/SCRI-main/supermassa_vcf2sm/ploidy_estimation/")
# res<-system("C:/Python27/python.exe ../supermassa/src/SuperMASSA.py --print_genotypes --inference f1 --ploidy_range 6:12 --file sugarcane_data/SugSNP225_progeny --f1_parent_data sugarcane_data/SugSNP225_parents", intern = TRUE)

# Auxiliary function
get_results <- function(res, progeny.file, parental.file){
  snp.name<-stringr::str_remove_all(progeny.file, "sugarcane_data/|_progeny")
  y0 <- res[stringr::str_which(res, "Prog_")]
  g <- t(sapply(stringr::str_split(y0, pattern = " \t\\(|,"), function(x) x[1:2]))
  y1 <- as.numeric(g[,2])
  names(y1) <- g[,1]
  x1 <- round(as.numeric(strsplit(res[3], split = " ")[[1]][8]), 5)
  x2 <- as.numeric(strsplit(res[4], split = " |,|\\(")[[1]][c(6,11,16)])
  z <- read.table(parental.file)[,1:3]
  names(x2) <- c("ploidy", sort(unique(z[,1])))
  y <- read.table(progeny.file, row.names = 1)
  z <- cbind(z, x2[z[,1]])
  df <- tibble::tibble(ind = rownames(y), 
                       A = y[,1],
                       a = y[,2],
                       dose = as.factor(y1[rownames(y)]))
  colnames(z) <- colnames(df)
  z <- z %>% group_by(ind) %>% summarise(across(everything(), list(mean)), .groups = 'drop')
  colnames(z) <- colnames(df)
  df <- rbind(df, z)
  p<-ggplot(df, aes(A, a)) + 
    xlim(min(df[,2:3]), max(df[,2:3])) +
    ylim(min(df[,2:3]), max(df[,2:3]))
  ploidy.lines<- c(c(0:x2[1]/x2[1]:0)[-(x2[1]+1)], 10e50)        
  for(l in 1:(x2[1]+1))  p<-p + geom_abline(intercept =0,  
                                            slope =ploidy.lines[l], 
                                            linetype = "dashed", 
                                            color = "gray")
  print(p + geom_point(aes(color = dose)) + 
          geom_point(aes(x = as.numeric(z[1,2]), y = as.numeric(z[1,3])), shape = 2, size = 2) +
          geom_point(aes(x = as.numeric(z[2,2]), y = as.numeric(z[2,3])), shape = 2, size = 2) +
          annotate(geom="text", x=max(df[,2:3]), 
                   y=max(df[,2:3]), label=x2[1], size = 7) +
          annotate(geom="text", x=max(df[,2:3]), 
                   y=max(df[,2:3]) - 3, label=round(x1,4), size = 5) +
          ggtitle(snp.name))
  return(data.frame(Pr = x1, ploidy = x2[1], "IACSP93-3046" = x2[2], "IACSP95-3018" = x2[3], row.names = snp.name))
} 

# Get results
get_results(res = res, 
            progeny.file = "sugarcane_data/SugSNP225_progeny", 
            parental.file = "sugarcane_data/SugSNP225_parents")

# Running in batch
fl <- list.files("~/repos/SCRI/supermassa_vcf2sm/ploidy_estimation/sugarcane_data")
fl

# Parental files
flpar <- grep(pattern = "par", fl, value = TRUE)
# Progeny files
flpro <- grep(pattern = "pro", fl, value = TRUE)

flpar
flpro

length(flpro)

pdf(file = "sugar_plots.pdf", width = 5, height = 5)
# Running in batch (~60 seconds)
df <- NULL
for(i in 1:length(flpro)){
  cat(i, "\n")
  src <- paste0("python2.7 ../supermassa/src/SuperMASSA.py --print_genotypes --inference f1 --ploidy_range 6:12 --file sugarcane_data/", 
                flpro[i], 
                " --f1_parent_data sugarcane_data/",
                flpar[i])
  res <- system(src, intern = TRUE)
  temp <- get_results(res, 
                    progeny.file = paste0("sugarcane_data/", flpro[i]), 
                    parental.file = paste0("sugarcane_data/",flpar[i]))
  df <- rbind(df, temp)
}
dev.off()
df
barplot(table(df[,2]), col = heat.colors(4))
