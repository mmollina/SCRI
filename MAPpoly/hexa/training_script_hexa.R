## Loading MAPpoly
require(mappoly)
setwd("~/repos/SCRI/MAPpoly/hexa/")

## Reading sweetpotato VCF files (chromosomes 3, 9 and 12)
dat <- NULL
for(i in c(3,9,12)){ # for all chromosomes use 1:15 (requires a lot of memory)
  cat("Loading chromosome", i, "...")
    tempfl <- tempfile(pattern = paste0("ch", i), fileext = ".vcf.gz")
    x <- "https://github.com/mmollina/MAPpoly_vignettes/raw/master/data/BT/sweetpotato_chr"
    address <- paste0(x, i, ".vcf.gz")
    download.file(url = address, destfile = tempfl, quiet = TRUE)
    dattemp <- read_vcf(file = tempfl, parent.1 = "PARENT1", parent.2 = "PARENT2", 
                        ploidy = 6, verbose = FALSE)
    dat <- merge_datasets(dat, dattemp)
  cat("\n")
}

## Plot original data
plot(dat)

## Filtering dataset by marker
dat <- filter_missing(input.data = dat, type = "marker", 
                      filter.thres = 0.05, inter = TRUE)
plot(dat)

## Filtering dataset by individual
dat <- filter_missing(input.data = dat, type = "individual", 
                      filter.thres = 0.05, inter = TRUE)

## Plot and print filtered data
plot(dat)
print(dat, detailed = TRUE)

## Segregation test
pval.bonf <- 0.05/dat$n.mrk ## approx Bonferroni correction
mrks.chi.filt <- filter_segregation(dat, 
                                    chisq.pval.thres =  pval.bonf, 
                                    inter = TRUE)
## Make initial sequence
seq.init<-make_seq_mappoly(mrks.chi.filt)
seq.init
plot(seq.init)
print(seq.init, detailed = TRUE)

## Two-point analysis
#~28 minutes using 24 CPUs (24 x Intel(R) Xeon(R) CPU E5-2670 v3 @ 2.30GHz)
tpt <- est_pairwise_rf(seq.init, ncpus = 24)
save(tpt, file = "~/repos/SCRI/MAPpoly/hexa/tpt_6602_mrk.rda")
## Two-point object to matrix
#~2.5 minutes
m <- rf_list_to_matrix(tpt, ncpus = 10)
save(m, file = "~/repos/SCRI/MAPpoly/hexa/m_6602_mrk.rda")

#load(file = "~/repos/SCRI/MAPpoly/hexa/tpt_6602_mrk.rda")
#load(file = "~/repos/SCRI/MAPpoly/hexa/m_6602_mrk.rda")

## Plot recombination fraction matrix
plot(m, fact = 10)

## Filtering markers by recombination fraction pairs
# Visualizing filtered matrix
mtemp <- rf_list_to_matrix(tpt, 
                           thresh.LOD.ph = 5, 
                           thresh.LOD.rf = 5,
                           thresh.rf = 0.15,
                           ncpus = 10)
plot(mtemp, fact = 10)
# Filtering
sf<-rf_snp_filter(tpt, 
                  thresh.LOD.ph = 5, 
                  thresh.LOD.rf = 5,
                  thresh.rf = 0.15, 
                  probs = c(0.05, 0.99))
mf<-make_mat_mappoly(m, sf)
plot(mf, fact = 10)

## Grouping
gr<-group_mappoly(mf, expected.groups = 3, comp.mat = TRUE)
gr

## Making ordered sequences using MDS for each group
## LG 1: corresponds to ch 3
s1 <- make_seq_mappoly(gr, 1, genomic.info = 1)
tpt1 <- make_pairs_mappoly(tpt, s1)
m1 <- rf_list_to_matrix(tpt1)
o1<-mds_mappoly(m1)
plot(o1)
so1 <- make_seq_mappoly(o1)
plot(m1, ord = so1$seq.mrk.names, fact = 10)
## Save data for ch 3
save(dat, o1, m1, so1, tpt1, file = "~/repos/SCRI/MAPpoly/hexa/ch3.rda")

## LG 2: corresponds to ch 9
s2 <- make_seq_mappoly(gr, 2, genomic.info = 1)
tpt2 <- make_pairs_mappoly(tpt, s2)
m2 <- rf_list_to_matrix(tpt2)
o2<-mds_mappoly(m2)
plot(o2)
so2 <- make_seq_mappoly(o2)
plot(m2, ord = so2$seq.mrk.names, fact = 10)
## Save data for ch 9
save(dat, o2, m2, so2, tpt2, file = "~/repos/SCRI/MAPpoly/hexa/ch9.rda")

## LG 3: corresponds to ch 12
s3 <- make_seq_mappoly(gr, 3, genomic.info = 1)
tpt3 <- make_pairs_mappoly(tpt, s3)
m3 <- rf_list_to_matrix(tpt3)
o3<-mds_mappoly(m3)
plot(o3)
so3 <- make_seq_mappoly(o3)
plot(m3, ord = so3$seq.mrk.names, fact = 10)
## Save data for ch 12
save(dat, o3, m3, so3, tpt3, file = "~/repos/SCRI/MAPpoly/hexa/ch12.rda")

## Run individual scripts for each linkage group
## R CMD BATCH --vanilla phase_ch3.R 
## R CMD BATCH --vanilla phase_ch9.R 
## R CMD BATCH --vanilla phase_ch12.R 
##
## After running the phasing scripts, open post_phasing.R and proceed 
## with the final steps



