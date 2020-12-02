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

X<-data.frame(snp_name = dat$mrk.names, 
              P1 = dat$dosage.p, 
              P2 = dat$dosage.q, 
              sequence = dat$sequence, 
              sequence_position = dat$sequence.pos,
              dat$geno.dose)
write.csv(X, file = "../data/B2721_dose.csv", row.names = FALSE)

Y <- read_geno_csv(file.in = "../data/B2721_dose.csv", ploidy = 4)
