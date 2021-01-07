setwd("~/repos/SCRI/supermassa_vcf2sm/")
fl <- list.files("sugarcane_data", full.names = TRUE)
flpar <- grep(pattern = "parents", fl, value = TRUE)
flpro <- grep(pattern = "progeny", fl, value = TRUE)

flpar[1]
flpro[1]

src1 <- "python supermassa/src/SuperMASSA.py --print_genotypes --inference f1 --ploidy_range 2:12"
src2 <- paste("--file", flpro[1])
src3 <- paste("--f1_parent_data", flpro[1])

system(paste(src1, src2, src3))
