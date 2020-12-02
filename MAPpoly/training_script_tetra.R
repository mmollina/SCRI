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
                    verbose = TRUE)
## Get genome information
source("get_solcap_snp_pos.R")
plot(dat)
print(dat, detailed = TRUE)

dev.off()
plot_mrk_info(dat, 738)
print_mrk(dat, 738)
plot_mrk_info(dat, "solcap_snp_c1_3722")
print_mrk(dat, "solcap_snp_c1_3722")

dat <- filter_missing(dat, type = "marker", 
                      filter.thres = 0.05)
dat <- filter_missing(input.data = dat, type = "individual", 
                      filter.thres = 0.025)
plot(dat)
s <- filter_segregation(input.data = dat, chisq.pval.thres = 0.05/dat$n.mrk)
s <- make_seq_mappoly(s)
plot(s)

#9.8 minutes in a Intel(R) Xeon(R) CPU E5-2670 v3 @ 2.30GHz (12 cores, 24 threads) and 128 GB RAM
nc <- parallel::detectCores() - 1
system.time(tpt <- est_pairwise_rf(s, ncpus = nc))
#tpt <- readRDS("tpt_6602_mrk.rds")
tpt$pairwise[1:10]
tpt$pairwise$`1-10`
plot(tpt, first.mrk = 1, second.mrk = 10)
#~1.1 minutes
system.time(m <- rf_list_to_matrix(tpt, ncpus = nc))
#m <- readRDS("m_6602_mrk.rds")
gen.ord <- get_genomic_order(s)
s.gen.ord <- make_seq_mappoly(gen.ord)
plot(m, ord = s.gen.ord$seq.mrk.names, fact = 10)
sf <- rf_snp_filter(input.twopt = tpt, 
                    thresh.LOD.ph = 5, 
                    thresh.LOD.rf = 5, 
                    thresh.rf = 0.25, 
                    probs = c(0.05, 0.95))
mf <- make_mat_mappoly(m, sf)
## Grouping
gr <- group_mappoly(mf, expected.groups = 12, 
                    comp.mat = TRUE)
gr

#### Functions ####
phasing_and_hmm_rf <- function(X){
  fl <- paste0("output_map_ch_", X$seq$sequence[1], ".txt")
  sink(fl)
  map <- est_rf_hmm_sequential(input.seq = X$seq,
                               start.set = 3,
                               thres.twopt = 10,
                               thres.hmm = 50,
                               extend.tail = 30,
                               twopt = X$tpt,
                               verbose = TRUE,
                               phase.number.limit = 20,
                               sub.map.size.diff.limit = 5) 
  sink()
  return(map)
}
error_model <- function(X, error = 0.05, tol = 10e-4){
  X$maps[[1]]$seq.rf <- rep(0.01,
                            length(X$maps[[1]]$seq.rf))
  x <- est_full_hmm_with_global_error(input.map = X, 
                                      error = error, 
                                      tol = tol, 
                                      verbose = FALSE)
  return(x)
}
#### Correspondence with genome
z<-as.numeric(colnames(gr$seq.vs.grouped.snp)[1:12])
#### Assembling linkage groups (order based on genome)
LGS<-vector("list", 12)
for(ch in 1:12){
  cat("\n ~~~~~~ ch:", ch, "...\n")
  lg <- which(z==ch)
  s.temp<-make_seq_mappoly(gr, lg, genomic.info = 1)
  tpt.temp <- make_pairs_mappoly(tpt, s.temp)
  s.temp.filt <- rf_snp_filter(tpt.temp, 5, 5, 0.15, c(0.05, 1))
  m.temp <- make_mat_mappoly(m, s.temp.filt)
  g.o <- get_genomic_order(s.temp)
  s.g <- make_seq_mappoly(g.o)
  tpt.temp <- make_pairs_mappoly(tpt, input.seq = s.g)
  LGS[[ch]] <- list(seq = s.g, tpt = tpt.temp)
}
#### Parallel map construction
cl <- parallel::makeCluster(12)
parallel::clusterEvalQ(cl, require(mappoly))
parallel::clusterExport(cl, "dat")
# ~12.5 minutes
MAPs.geno <- parallel::parLapply(cl, LGS, phasing_and_hmm_rf)
# ~2.5 minutes
MAPs.geno.err <- parallel::parLapply(cl, MAPs.geno, 
                                     error_model)
# ~22 seconds
genoprob <- parallel::parLapply(cl,
                                MAPs.geno.err,
                                calc_genoprob_error, 
                                step = 1, 
                                error = 0.05)
parallel::stopCluster(cl)

#### Map results
map.out <- plot_map_list(MAPs.geno.err, col = "ggstyle")
summary_maps(MAPs.geno.err)

## trying other arguments for LG 4
temp.map <- MAPs.geno.err[[4]]
plot(temp.map)
s.temp <- make_seq_mappoly(temp.map)
tpt.temp <- make_pairs_mappoly(tpt, s.temp)
m.temp <- rf_list_to_matrix(tpt.temp)
plot(m.temp, ord = s.temp$seq.mrk.names, fact = 1)
lg4.geno.map<-est_rf_hmm_sequential(input.seq = s.temp,
                                    start.set = 3,
                                    thres.twopt = 10,
                                    thres.hmm = 50,
                                    extend.tail = 100,
                                    twopt = tpt.temp,
                                    verbose = TRUE,
                                    tol = 10e-2,
                                    tol.final = 10e-4,
                                    phase.number.limit = 20,
                                    sub.map.size.diff.limit = 2,
                                    info.tail = TRUE)
lg4.geno.map.err <- est_full_hmm_with_global_error(lg4.geno.map, 
                                                   error = 0.05, 
                                                   tol = 10e-4, 
                                                   verbose = TRUE)
plot(lg4.geno.map.err)
#MAPs.geno.err[[4]] <- lg4.geno.map

## dropping initial markers in LG 12
temp.map <- MAPs.geno.err[[12]]
plot(temp.map)
s.temp <- make_seq_mappoly(temp.map)
tpt.temp <- make_pairs_mappoly(tpt, s.temp)
m.temp <- rf_list_to_matrix(tpt.temp)
plot(m.temp, ord = s.temp$seq.mrk.names)
x <- split_and_rephase(temp.map, twopt = tpt.temp, gap.threshold = 8)
plot(x)
plot(temp.map, left.lim = 0, right.lim = 20, mrk.names = TRUE)
temp.map <- drop_marker(temp.map, 
                        mrk = c("solcap_snp_c2_24586", 
                                "solcap_snp_c2_24565", 
                                "solcap_snp_c2_24560"))
lg12.geno.map <- est_full_hmm_with_global_error(temp.map, 
                                                error = 0.05, 
                                                tol = 10e-4, 
                                                verbose = TRUE)
plot(lg12.geno.map)

MAPs.geno.err[[12]] <- lg12.geno.map


## Final map
final.map <- update_map(MAPs.geno.err)
map.out <- plot_map_list(MAPs.geno.err, col = "ggstyle")
map.out
plot_genome_vs_map(final.map, same.ch.lg = TRUE)
summary_maps(final.map)
export_map_list(final.map, file = "output_map.csv")

#### Preferential pairing
pp <- calc_prefpair_profiles(genoprob)
print(pp)
head(pp$prefpair.psi)
plot(pp, P = "Atlantic", Q = "B1829")

#### Haplotype probabilities
hp <- calc_homoprob(genoprob)
print(hp)
plot(hp, ind = 2, lg = 1)
plot(hp, ind = 2, lg = 1:12, use.plotly = FALSE)



plot(final.map[[1]])
plot(final.map[[2]], xlim = c(0, 154.07))
plot(final.map[[3]], xlim = c(0, 154.07))
plot(final.map[[4]], xlim = c(0, 154.07))
plot(final.map[[5]], xlim = c(0, 154.07))
plot(final.map[[6]], xlim = c(0, 154.07))
plot(final.map[[7]], xlim = c(0, 154.07))
plot(final.map[[8]], xlim = c(0, 154.07))
plot(final.map[[9]], xlim = c(0, 154.07))
plot(final.map[[10]], xlim = c(0, 154.07))
plot(final.map[[11]], xlim = c(0, 154.07))
plot(final.map[[12]], xlim = c(0, 154.07))



save.image(file = "image.rda")









#### Functions ####
phasing_and_hmm_rf <- function(X){
  fl <- paste0("output_map_ch_", X$seq$sequence[1], ".txt")
  sink(fl)
  map <- est_rf_hmm_sequential(input.seq = X$seq,
                               start.set = 3,
                               thres.twopt = 10,
                               thres.hmm = 50,
                               extend.tail = 30,
                               twopt = X$tpt,
                               verbose = TRUE,
                               phase.number.limit = 20,
                               sub.map.size.diff.limit = 5) 
  sink()
  return(map)
}
error_model <- function(X, error = 0.05, tol = 10e-4){
  X$maps[[1]]$seq.rf <- rep(0.01,
                            length(X$maps[[1]]$seq.rf))
  x <- est_full_hmm_with_global_error(input.map = X, 
                                      error = error, 
                                      tol = tol, 
                                      verbose = FALSE)
  return(x)
}
#### Correspondence with genome
z<-as.numeric(colnames(gr$seq.vs.grouped.snp)[1:12])
#### Assembling linkage groups (order based on MDS)
MDS.o <- LGS.mds<-vector("list", 12)
#op<-par(ask = TRUE)
for(ch in 1:12){
  cat("\n ~~~~~~ ch:", ch, "...\n")
  lg <- which(z==ch)
  s.temp<-make_seq_mappoly(gr, lg, genomic.info = 1)
  tpt.temp <- make_pairs_mappoly(tpt, s.temp)
  s.temp.filt <- rf_snp_filter(tpt.temp, 5, 5, 0.15, c(0.05, 1))
  m.temp <- make_mat_mappoly(m, s.temp.filt)
  MDS.o[[ch]] <- mds_mappoly(m.temp)
  s.mds <- make_seq_mappoly(MDS.o[[ch]])
  plot(m.temp, ord = s.mds$seq.mrk.names)
  tpt.temp <- make_pairs_mappoly(tpt, input.seq = s.mds)
  LGS.mds[[ch]] <- list(seq = s.mds, tpt = tpt.temp)
}
#par(op)
#### Parallel map construction
cl <- parallel::makeCluster(12)
parallel::clusterEvalQ(cl, require(mappoly))
parallel::clusterExport(cl, "dat")
# ~32.2 minutes
system.time(MAPs.mds <- parallel::parLapply(cl, LGS.mds, phasing_and_hmm_rf))
# ~3.2 minutes
system.time(MAPs.mds.err <- parallel::parLapply(cl, MAPs.mds, 
                                     error_model))
# ~17 seconds
system.time(genoprob.mds <- parallel::parLapply(cl,
                                MAPs.mds.err,
                                calc_genoprob_error, 
                                step = 1, 
                                error = 0.05))
parallel::stopCluster(cl)
plot_map_list(MAPs.mds.err)
summary_maps(MAPs.mds.err)
plot_genome_vs_map(MAPs.mds.err, same.ch.lg = TRUE)


## dropping initial markers in LG 5
temp.map <- MAPs.mds.err[[5]]
plot(temp.map)
s.temp <- make_seq_mappoly(temp.map)
tpt.temp <- make_pairs_mappoly(tpt, s.temp)
m.temp <- rf_list_to_matrix(tpt.temp)
plot(m.temp, ord = s.temp$seq.mrk.names)
plot(temp.map, left.lim = 0, right.lim = 80, mrk.names = TRUE)
plot(temp.map, left.lim = 30, right.lim = 78, mrk.names = TRUE)
plot(temp.map, left.lim = 35, right.lim = 50, mrk.names = TRUE)
temp.map <- drop_marker(temp.map, mrk = c("solcap_snp_c1_14797", "solcap_snp_c2_50313"))
plot(temp.map)
lg5.geno.map <- est_full_hmm_with_global_error(temp.map, 
                                                error = 0.05, 
                                                tol = 10e-4, 
                                                verbose = TRUE)
plot(lg5.geno.map)
x <- split_and_rephase(lg5.geno.map, twopt = tpt.temp, gap.threshold = 5)
plot(x)
lg5.geno.map <- est_full_hmm_with_global_error(x, 
                                               error = 0.05, 
                                               tol = 10e-4, 
                                               verbose = TRUE)
plot(lg5.geno.map)
MAPs.mds.err[[5]] <- lg5.geno.map
plot_map_list(MAPs.mds.err)
plot_map_list(final.map)
summary_maps(MAPs.mds.err)
plot_genome_vs_map(MAPs.mds.err, same.ch.lg = TRUE)

map0 <- MAPs.mds.err[[5]]
m.temp <- make_mat_mappoly(m, make_seq_mappoly(map0))
map1 <- reest_rf(MAPs.mds.err[[5]], input.mat = m.temp, method = "ols")
map2 <- reest_rf(MAPs.mds.err[[5]], input.mat = m.temp, method = "wMDS_to_1D_pc", input.mds = MDS.o[[5]])

plot_map_list(list(map0, map1, map2))
