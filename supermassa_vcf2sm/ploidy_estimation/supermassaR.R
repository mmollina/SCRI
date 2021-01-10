# Using 'system'
setwd("~/repos/SCRI/supermassa_vcf2sm/ploidy_estimation/")
system("python ../supermassa/src/SuperMASSA.py --print_genotypes --inference f1 --ploidy_range 2:12 --naive_posterior_reporting_threshold 0.8 --file SugSNP_382.pro --f1_parent_data SugSNP_382.par", 
       intern = TRUE)
# Running in batch
setwd("~/repos/SCRI/supermassa_vcf2sm/ploidy_estimation/sugarcane_data/")
fl <- list.files("~/repos/SCRI/supermassa_vcf2sm/ploidy_estimation/sugarcane_data")
fl
flpar <- grep(pattern = "par", fl, value = TRUE)
flpro <- grep(pattern = "pro", fl, value = TRUE)

flpar[1]
flpro[1]

df <- NULL
for(i in 1:length(flpro)){
  cat(flpro[i], "\n")
  src1 <- "python ../../supermassa/src/SuperMASSA.py --print_genotypes --inference f1 --ploidy_range 2:12 --naive_posterior_reporting_threshold 0.8"
  src2 <- paste0("--file ~/repos/SCRI/supermassa_vcf2sm/ploidy_estimation/sugarcane_data/", flpro[i])
  src3 <- paste0("--f1_parent_data ~/repos/SCRI/supermassa_vcf2sm/ploidy_estimation/sugarcane_data/", flpro[i])
  src <- paste(src1, src2, src3)
  res <- system(src, intern = TRUE)
  x1 <- round(as.numeric(strsplit(res[3], split = " ")[[1]][8]), 5)
  x2 <- as.numeric(strsplit(res[4], split = " |,|\\(")[[1]][c(6,11,16)])
  df <- rbind(df, c(x1,x2))
}
colnames(df) <- c("post.ploidy", "ploidy", "d.P1", "d.P2")
rownames(df) <- stringr::str_remove_all(flpar, "_parents")
df
barplot(table(df[,2]), col = heat.colors(4))
