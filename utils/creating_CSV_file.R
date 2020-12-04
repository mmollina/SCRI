require(mappoly)
setwd("~/repos/SCRI/MAPpoly/")
## Downloading and reading B2721 fitpolyprobabilistic scores
address <- "https://github.com/mmollina/SCRI/raw/main/data/fitpoly_tetra_call/B2721_scores.zip"
tempfl <- tempfile(pattern = "B2721_CC_", fileext = ".zip")
download.file(url = address, destfile = tempfl)
unzip(tempfl, files = "B2721_scores.dat")
unlink(tempfl)
##96 seconds
dat <- read_fitpoly(file.in = "B2721_scores.dat", 
                    ploidy = 4, 
                    parent1 = "Atlantic", 
                    parent2 = "B1829", 
                    verbose = TRUE, 
                    filter.non.conforming = FALSE, 
                    elim.redundant = FALSE, 
                    prob.thres = 0.9)
source("get_solcap_snp_pos.R")
X<-data.frame(snp_name = dat$mrk.names, 
              P1 = dat$dosage.p, 
              P2 = dat$dosage.q, 
              sequence = dat$sequence, 
              sequence_position = dat$sequence.pos,
              dat$geno.dose)
write.csv(X, file = "../data/B2721_dose.csv", row.names = FALSE)
Y <- read_geno_csv(file.in = "../data/B2721_dose.csv", ploidy = 4)


dat <- read_fitpoly(file.in = "B2721_scores.dat", 
                    ploidy = 4, 
                    parent1 = "Atlantic", 
                    parent2 = "B1829", 
                    verbose = TRUE,
                    prob.thres = 0.9)

## Writing mappoly input file
setwd("../data/")
indnames<-as.character(dat$ind.names)
mrknames<-as.character(dat$mrk.names)
write(paste("ploidy", dat$m), file="B2721_mappoly_dose")
write(paste("nind", length(indnames)), file="B2721_mappoly_dose", append=TRUE)
write(paste("nmrk", length(mrknames)), file="B2721_mappoly_dose", append=TRUE)
cat("mrknames", mrknames, file="B2721_mappoly_dose", append=TRUE)
cat("\nindnames", indnames, file="B2721_mappoly_dose", append=TRUE)
cat("\ndosageP", dat$dosage.p[mrknames], file="B2721_mappoly_dose", append=TRUE)
cat("\ndosageQ", dat$dosage.q[mrknames], file="B2721_mappoly_dose", append=TRUE)
cat("\nseq", dat$sequence, file="B2721_mappoly_dose", append=TRUE)
cat("\nseqpos", dat$sequence.pos, file="B2721_mappoly_dose", append=TRUE)
write("\nnphen 0", file="B2721_mappoly_dose", append=TRUE)
write("pheno---------------------------------------", file="B2721_mappoly_dose", append=TRUE)
write("geno---------------------------------------", file="B2721_mappoly_dose", append=TRUE)
write.table(dat$geno.dose, file="B2721_mappoly_dose", append=TRUE, quote=FALSE,
            row.names=FALSE, col.names=FALSE)
dat.dose.mpl <- read_geno("B2721_mappoly_dose")
save(dat.dose.mpl, file = "B2721_mappoly_dose.RData", compress = "xz", compression_level = 9)


## Writing mappoly input file
setwd("../data/")
indnames<-as.character(dat$ind.names)
mrknames<-as.character(dat$mrk.names)
write(paste("ploidy", dat$m), file="B2721_mappoly_prob")
write(paste("nind", length(indnames)), file="B2721_mappoly_prob", append=TRUE)
write(paste("nmrk", length(mrknames)), file="B2721_mappoly_prob", append=TRUE)
cat("mrknames", mrknames, file="B2721_mappoly_prob", append=TRUE)
cat("\nindnames", indnames, file="B2721_mappoly_prob", append=TRUE)
cat("\ndosageP", dat$dosage.p[mrknames], file="B2721_mappoly_prob", append=TRUE)
cat("\ndosageQ", dat$dosage.q[mrknames], file="B2721_mappoly_prob", append=TRUE)
cat("\nseq", dat$sequence, file="B2721_mappoly_prob", append=TRUE)
cat("\nseqpos", dat$sequence.pos, file="B2721_mappoly_prob", append=TRUE)
write("\nnphen 0", file="B2721_mappoly_prob", append=TRUE)
write("pheno---------------------------------------", file="B2721_mappoly_prob", append=TRUE)
write("geno---------------------------------------", file="B2721_mappoly_prob", append=TRUE)
write.table(dat$geno, file="B2721_mappoly_prob", append=TRUE, quote=FALSE,
            row.names=FALSE, col.names=FALSE)
dat.prob.mpl <- read_geno_prob("B2721_mappoly_prob")
save(dat.prob.mpl, file = "B2721_mappoly_prob.RData", compress = "xz", compression_level = 9)
