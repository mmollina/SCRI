require(mappoly)
dat <- NULL
for(i in c(3,9,12)){
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
plot(dat)
## Filtering dataset by marker
dat <- filter_missing(input.data = dat, type = "marker", 
                      filter.thres = 0.05, inter = TRUE)
plot(dat)
## Filtering dataset by individual
dat <- filter_missing(input.data = dat, type = "individual", 
                      filter.thres = 0.05, inter = TRUE)
plot(dat)
print(dat, detailed = TRUE)

## Segregation test
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, 
                                    chisq.pval.thres =  pval.bonf, 
                                    inter = TRUE)
seq.init<-make_seq_mappoly(mrks.chi.filt)
seq.init
plot(seq.init)
print(seq.init, detailed = TRUE)

#~28 minutes
tpt <- est_pairwise_rf(seq.init, ncpus = 24)
saveRDS(tpt, file = "~/repos/SCRI_training/MAPpoly/tpt_6602_mrk.rds")
#~2.5 minutes
m <- rf_list_to_matrix(tpt, ncpus = 10)
saveRDS(m, file = "~/repos/SCRI_training/MAPpoly/m_6602_mrk.rds")

tpt <- readRDS(file = "~/repos/SCRI_training/MAPpoly/tpt_6602_mrk.rds")
m <- readRDS(file = "~/repos/SCRI_training/MAPpoly/m_6602_mrk.rds")

plot(m, fact = 10)

mtemp <- rf_list_to_matrix(tpt, 
                        thresh.LOD.ph = 5, 
                        thresh.LOD.rf = 5,
                        thresh.rf = 0.15,
                        ncpus = 10)
plot(mtemp, fact = 10)

sf<-rf_snp_filter(tpt, 
                  thresh.LOD.ph = 5, 
                  thresh.LOD.rf = 5,
                  thresh.rf = 0.15, 
                  probs = c(0.05, 0.99))

mf<-make_mat_mappoly(m, sf)
plot(mf, fact = 10)
gr<-group_mappoly(mf, expected.groups = 3, comp.mat = TRUE)

s1<-make_seq_mappoly(gr, 1, genomic.info = 1)
tpt1<-make_pairs_mappoly(tpt, s1)
st<-rf_snp_filter(tpt1, 5, 5, 0.15, c(0.1, 1))
m1 <- make_mat_mappoly(input.mat = mf, st)
plot(m1, ord = st$seq.mrk.names, fact = 10)
o1<-mds_mappoly(m1, n = c(1685, 733, 786, 383))
plot(o1)
so <- make_seq_mappoly(o1)
plot(m1, ord = so$seq.mrk.names, fact = 10)


## Phasing and re-estimating genetic map
system.time(lg1.map<-est_rf_hmm_sequential(input.seq = st,
                                           start.set = 3,
                                           thres.twopt = 10,
                                           thres.hmm = 20,
                                           extend.tail = 50,
                                           twopt = tpt1,
                                           verbose = TRUE,
                                           tol = 10e-2,
                                           tol.final = 10e-3,
                                           phase.number.limit = 20,
                                           sub.map.size.diff.limit =  5,
                                           info.tail = TRUE,
                                           reestimate.single.ph.configuration = TRUE))

plot(so$sequence.pos)












require(mappoly)
dat <- NULL
for(i in c(3,9,12)){
  cat("Loading chromosome", i, "...\n")
  invisible(capture.output(y <- {
    tempfl <- tempfile(pattern = paste0("ch", i), fileext = ".vcf.gz")
    x <- "https://github.com/mmollina/MAPpoly_vignettes/raw/master/data/BT/sweetpotato_chr"
    address <- paste0(x, i, ".vcf.gz")
    download.file(url = address, destfile = tempfl)
    dattemp <- read_vcf(file = tempfl, parent.1 = "PARENT1", parent.2 = "PARENT2", ploidy = 6)
    dat <- merge_datasets(dat, dattemp)
  }))
  cat("\n")
}
## Filtering dataset by marker
dat <- filter_missing(input.data = dat, type = "marker", 
                      filter.thres = 0.05, inter = TRUE)

## Filtering dataset by individual
dat <- filter_missing(input.data = dat, type = "individual", 
                      filter.thres = 0.05, inter = TRUE)
print(dat, detailed = TRUE)

## Segregation test
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, 
                                    chisq.pval.thres =  pval.bonf, 
                                    inter = TRUE)
seq.init<-make_seq_mappoly(mrks.chi.filt)
length(seq.init$seq.mrk.names)
plot(seq.init)
print(seq.init, detailed = TRUE)

z<-system.time(tpt <- est_pairwise_rf(seq.init, ncpus = 14))
z[3]/3600
m<-rf_list_to_matrix(tpt)
sf<-rf_snp_filter(tpt, 5,5, 0.15, probs = c(0.05, 0.99))
mf<-make_mat_mappoly(m, sf)
plot(mf, fact = 5)
gr<-group_mappoly(mf, expected.groups = 3, comp.mat = TRUE)
s1<-make_seq_mappoly(gr, 2, genomic.info = 1)
tpt1<-make_pairs_mappoly(tpt, s1)
st<-rf_snp_filter(tpt1, 10, 10, 0.25, c(0.1, 1))
m1 <- make_mat_mappoly(input.mat = mf, st)
plot(m1, ord = st$seq.mrk.names, fact = 5)
#o1<-mds_mappoly(m1)
#plot(o1)
#so <- make_seq_mappoly(o1)
#plot(m1, ord = so$seq.mrk.names, fact = 5)
## Phasing and re-estimating genetic map
system.time(lg1.map<-est_rf_hmm_sequential(input.seq = st,
                                           start.set = 3,
                                           thres.twopt = 10,
                                           thres.hmm = 20,
                                           extend.tail = 50,
                                           twopt = tpt1,
                                           verbose = TRUE,
                                           tol = 10e-2,
                                           tol.final = 10e-3,
                                           phase.number.limit = 20,
                                           sub.map.size.diff.limit =  5,
                                           info.tail = TRUE,
                                           reestimate.single.ph.configuration = TRUE))

plot(so$sequence.pos)





