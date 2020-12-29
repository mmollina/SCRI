load("~/repos/SCRI/MAPpoly/tetra/all_tetra_analysis.rda")
library("mappoly")
setwd("~/repos/SCRI/MAPpoly/tetra/")

#### Downloading data from GitHub (~9 seconds) #### 
address <- "https://github.com/mmollina/SCRI/raw/main/data/fitpoly_tetra_call/B2721_scores.zip"
tempfl <- tempfile(pattern = "B2721_", fileext = ".zip")
download.file(url = address, destfile = tempfl)
unzip(tempfl, files = "B2721_scores.dat")

#### Reading fitPoly data (~100 seconds) #### 
dat <- read_fitpoly(file.in = "B2721_scores.dat", 
                    ploidy = 4, 
                    parent1 = "Atlantic",
                    parent2 = "B1829", 
                    prob.thres = 0.95,
                    verbose = TRUE)
class(dat)
plot(dat)
print(dat, detailed = TRUE)
## Adding genome information
source("get_solcap_snp_pos.R")
print(dat, detailed = TRUE)

## Inspecting a specific SNP
plot_mrk_info(dat, 738)
print_mrk(dat, 738)
segreg_poly(m = 4, dP = 1, dQ = 2)
#### Filtering ####
## Filtering by marker
dat <- filter_missing(input.data = dat, 
                      type = "marker", 
                      filter.thres = 0.05)

## Filter by individual
dat <- filter_missing(input.data = dat,
                      type = "individual", 
                      filter.thres = 0.025)
plot(dat)

## Filter by segregation
s <- filter_segregation(input.data = dat, chisq.pval.thres = 0.05/dat$n.mrk)
class(s)

#### Making a marker sequence ####
s <- make_seq_mappoly(s)
s
class(s)
plot(s)

s1_30 <- make_seq_mappoly(dat, 1:30)
s1_30
plot(s1_30)

s10 <- make_seq_mappoly(dat, "seq10")
s10
plot(s10)

#### Two-point analysis #### ~8 seconds
nc <- parallel::detectCores() - 1
nc
tpt10 <- est_pairwise_rf(s10, ncpus = nc, est.type = "disc")
tpt10
tpt10$pairwise[90:92]
round(tpt10$pairwise[[91]], 2)
plot(tpt10, first.mrk = 1, second.mrk = 2074)

# ~9.8 minutes
# tpt <- est_pairwise_rf(s, ncpus = nc)
# save(tpt, file = "two_point_rfs.rda")
# ~40 seconds
load("two_point_rfs.rda")

## Recombination Fraction Matrix
# ~ 80 seconds
# m <- rf_list_to_matrix(tpt, ncpus = nc)
# save(m, file = "rfs_mat.rda")
load("rfs_mat.rda")

# Plot recombination fraction matrix
gen.ord <- get_genomic_order(s)
gen.ord
plot(gen.ord)
s.gen.ord <- make_seq_mappoly(gen.ord)
plot(m, ord = s.gen.ord$seq.mrk.names, fact = 10)

#### Grouping ####
gr <- group_mappoly(m, 
                    expected.groups = 12, 
                    comp.mat = TRUE)
gr

## Make sequence from group
s10 <- make_seq_mappoly(gr, 
                       arg = 1, 
                       genomic.info = 1)
s10
tpt10 <- make_pairs_mappoly(tpt, s10)
m10 <- rf_list_to_matrix(tpt10)
plot(m10)

#### Ordering markers ####
## Genome order
geno.o.10 <- get_genomic_order(s10)
s.geno.o.10 <- make_seq_mappoly(geno.o.10)
plot(m10, ord = s.geno.o.10$seq.mrk.names)
## MDS order
mds.o.10 <- mds_mappoly(m10)
plot(mds.o.10)
s.mds.o.10 <- make_seq_mappoly(mds.o.10)
plot(m10, ord = s.mds.o.10$seq.mrk.names)

## Comparing genome vs MDS
plot(match(s.geno.o.10$seq.mrk.names, s.mds.o.10$seq.mrk.names), 
     col = 2, pch = 20, ylab = "MDS", xlab = "genome")


#### Filtering using recombination fractions ####
mf10 <- rf_list_to_matrix(tpt10, 
                          thresh.LOD.ph = 5, 
                          thresh.LOD.rf = 5,
                          thresh.rf = 0.15)
plot(mf10, ord = s.geno.o.10$seq.mrk.names)
## Filtering
s.10.f <- rf_snp_filter(input.twopt = tpt10, 
                    thresh.LOD.ph = 5, 
                    thresh.LOD.rf = 5, 
                    thresh.rf = 0.15, 
                    probs = c(0.05, 1))

s.geno.o.10
s.mds.o.10

#### Phasing using genome order ~ 2 minutes ####
lg10.geno.map <- est_rf_hmm_sequential(input.seq = s.geno.o.10,
                                       start.set = 3,
                                       thres.twopt = 10,
                                       thres.hmm = 50,
                                       extend.tail = 30,
                                       twopt = tpt10,
                                       verbose = TRUE,
                                       tol = 10e-2,
                                       tol.final = 10e-4,
                                       phase.number.limit = 20,
                                       sub.map.size.diff.limit = 5,
                                       info.tail = TRUE)
print(lg10.geno.map, detailed = TRUE)
plot(lg10.geno.map)
plot(lg10.geno.map, left.lim = 38, right.lim = 42, mrk.names = TRUE)

#### Phasing using MDS order ~2 minutes####
system.time(lg10.mds.map <- est_rf_hmm_sequential(input.seq = s.mds.o.10,
                                       start.set = 3,
                                       thres.twopt = 10,
                                       thres.hmm = 50,
                                       extend.tail = 30,
                                       twopt = tpt10,
                                       verbose = TRUE,
                                       tol = 10e-2,
                                       tol.final = 10e-4,
                                       phase.number.limit = 20,
                                       sub.map.size.diff.limit = 5,
                                       info.tail = TRUE))
print(lg10.mds.map, detailed = TRUE)
plot(lg10.mds.map)
plot_map_list(list(geno = lg10.geno.map, mds = lg10.mds.map), 
              col = c("turquoise3", "orangered2"))

#### Modeling errors ####
## ~60 seconds
lg10.geno.map.err<-est_full_hmm_with_global_error(lg10.geno.map, 
                                                  error = 0.05, 
                                                  verbose = TRUE)
## ~50 seconds
lg10.mds.map.err<-est_full_hmm_with_global_error(lg10.mds.map, 
                                                 error = 0.05, 
                                                 verbose = TRUE)
plot_map_list(list(geno = lg10.geno.map, 
                   geno.err = lg10.geno.map.err, 
                   mds = lg10.mds.map,
                   mds.err = lg10.mds.map.err), 
              col = c("turquoise4", "turquoise2", 
                      "orangered3", "orangered1"))

#### Genotype and homolog probabilities ####
## No error
g.lg10 <- calc_genoprob(lg10.geno.map, step = 1)
## Homolog probabilities
h.lg10 <- calc_homoprob(g.lg10)
plot(h.lg10, ind = "B2721.110")

## Error
g.lg10.err <- calc_genoprob_error(lg10.geno.map.err, 
                                  step = 1, 
                                  error = 0.05)
## Homolog probabilities
h.lg10.err <- calc_homoprob(g.lg10.err)
plot(h.lg10.err, ind = "B2721.110")

#### Preferential probability profiles ####
pp.lg10.err <- calc_prefpair_profiles(g.lg10.err)
plot(pp.lg10.err)

#### Parallel Mapping ####
#### Functions 
phasing_and_hmm_rf <- function(X){
  dir.create("map_output", showWarnings = FALSE)
  fl <- paste0("output_map_ch_", X$seq$sequence[1], ".txt")
  fl <- file.path("map_output", fl)
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
  m.temp <- make_mat_mappoly(m, s.temp)
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
plot_map_list(MAPs.geno, col = "ggstyle")
# ~2.5 minutes
MAPs.geno.err <- parallel::parLapply(cl, MAPs.geno, 
                                     error_model)
plot_map_list(MAPs.geno.err, col = "ggstyle")
# ~22 seconds
genoprob <- parallel::parLapply(cl,
                                MAPs.geno.err,
                                calc_genoprob_error, 
                                step = 1, 
                                error = 0.05)
plot_map_list(MAPs.geno.err)
parallel::stopCluster(cl)
#### Homolog and preferential pairing probabilities ####
## Homologs
hp <- calc_homoprob(genoprob)
print(hp)
plot(hp, ind = 5, lg = 3)
plot(hp, ind = 5, lg = 1:12, use.plotly = FALSE)
## Preferential pairing
pp <- calc_prefpair_profiles(genoprob)
print(pp)
head(pp$prefpair.psi)
plot(pp, P = "Atlantic", Q = "B1829")

#### Summary and exporting ####
final.map <- update_map(MAPs.geno.err)
summary_maps(final.map)
plot_genome_vs_map(final.map, 
                   same.ch.lg = TRUE)
plot_genome_vs_map(final.map, 
                   same.ch.lg = FALSE)
plot_map_list(final.map, horiz = FALSE, 
              col = viridis::viridis(12))
export_map_list(final.map, file = "map_output.csv")




