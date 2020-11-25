require(mappoly)

## Downloading and reading B2721 fitpolyprobabilistic scores
setwd("~/repos/SCRI_MAPpoly_training/")
address <- "https://github.com/mmollina/SCRI_MAPpoly_training/raw/main/data/B2721_scores.zip"
tempfl <- tempfile(pattern = "B2721_CC_", fileext = ".zip")
download.file(url = address, destfile = tempfl)
unzip(tempfl, files = "B2721_scores.dat")
B2721 <- read_fitpoly(file.in = "B2721_scores.dat", 
                      ploidy = 4, 
                      parent1 = "Atlantic", 
                      parent2 = "B1829", 
                      verbose = TRUE)
unlink(tempfl)

## Including genome position (using S. tuberosum genome v4.03)
read.csv("https://github.com/mmollina/SCRI_MAPpoly_training/raw/main/data/solcap_snp_pos_V4.03.csv")
B2721$sequence <- solcap.snp.pos[B2721$mrk.names,"chr"]
B2721$sequence.pos <- solcap.snp.pos[B2721$mrk.names,"pos"]
names(B2721$sequence.pos) <- names(B2721$sequence) <- B2721$mrk.names


plot(B2721)
print(B2721, detailed = TRUE)
plot_mrk_info(B2721, 37)
print_mrk(B2721, 37)
dat <- filter_missing(B2721, type = "marker", filter.thres = 0.05, inter = F)
dat <- filter_missing(dat, type = "individual", filter.thres = 0.025, inter = F)
plot(dat)
s <- filter_segregation(dat, chisq.pval.thres = 0.05/dat$n.mrk, inter = F)
s <- make_seq_mappoly(s)
plot(s)

#9.8 minutes in a Intel(R) Xeon(R) CPU E5-2670 v3 @ 2.30GHz (12 cores, 24 threads) and 128 GB RAM
nc <- parallel::detectCores() - 1
system.time(tpt <- est_pairwise_rf(s, ncpus = nc))
plot_mrk_info(B2721, 1)
plot_mrk_info(B2721, "solcap_snp_c1_1")
plot_mrk_info(B2721, 10)
plot_mrk_info(B2721, "solcap_snp_c1_10042")
tpt$pairwise[1:10]
tpt$pairwise$`1-10`
plot(tpt, first.mrk = 1, second.mrk = 10)
#~1.1 minutes
m <- rf_list_to_matrix(tpt, ncpus = 24)
plot(m, ord = rownames(get_genomic_order(s)), fact = 10)
mf <- rf_list_to_matrix(tpt, 
                        thresh.LOD.ph = 5, 
                        thresh.LOD.rf = 5, 
                        thresh.rf = 0.25)
plot(mf, ord = rownames(get_genomic_order(s)), fact = 10)
sf <- rf_snp_filter(input.twopt = tpt, 
                    thresh.LOD.ph = 5, 
                    thresh.LOD.rf = 5, 
                    thresh.rf = 0.25, 
                    probs = c(0.05, 0.95))

mf <- make_mat_mappoly(m, sf)

gr <- group_mappoly(mf, expected.groups = 12, comp.mat = TRUE)
gr


a <- match(1:12, as.numeric(colnames(gr$seq.vs.grouped.snp)[-13]))
L<-vector("list", 12)
for(i in 1:11){
  s1 <- make_seq_mappoly(gr, arg = a[i], genomic.info = 1)
  tpt1 <- make_pairs_mappoly(tpt, s1)
  m1 <- make_mat_mappoly(m, s1)
  ##Genomic order
  g.o1 <- get_genomic_order(s1)
  s1.g <- make_seq_mappoly(g.o1)
  plot(g.o1)
  plot(m1, ord = s1.g$seq.mrk.names, fact = 3)
  ##3.8 minutes
  system.time(lg1.geno.map<-est_rf_hmm_sequential(input.seq = s1.g,
                                                  start.set = 3,
                                                  thres.twopt = 10,
                                                  thres.hmm = 50,
                                                  extend.tail = 100, #30 for all, except 11
                                                  twopt = tpt1,
                                                  verbose = TRUE,
                                                  tol = 10e-2,
                                                  tol.final = 10e-4,
                                                  phase.number.limit = 20,
                                                  sub.map.size.diff.limit = 5,
                                                  info.tail = TRUE,
                                                  reestimate.single.ph.configuration = TRUE))
  plot(lg1.geno.map)
  ##3.3 minutes
  system.time(lg1.geno.map.err <- est_full_hmm_with_global_error(lg1.geno.map, error = 0.05, tol = 10e-4))
  plot(lg1.geno.map.err)
  ##MDS order
  mds.o1 <- mds_mappoly(m1)
  plot(mds.o1)
  s1.mds <- make_seq_mappoly(mds.o1)
  plot(m1, ord = s1.g$seq.mrk.names, fact = 3)
  ## 5.1 minutes
  system.time(lg1.mds.map<-est_rf_hmm_sequential(input.seq = s1.mds,
                                                 start.set = 3,
                                                 thres.twopt = 10,
                                                 thres.hmm = 50,
                                                 extend.tail = 30,
                                                 twopt = tpt1,
                                                 verbose = TRUE,
                                                 tol = 10e-2,
                                                 tol.final = 10e-4,
                                                 phase.number.limit = 20,
                                                 sub.map.size.diff.limit = 5,
                                                 info.tail = TRUE,
                                                 reestimate.single.ph.configuration = TRUE))
  plot(lg1.mds.map)
  ## 2.1 minutes
  system.time(lg1.mds.map.err <- est_full_hmm_with_global_error(lg1.mds.map, error = 0.05, tol = 10e-4))
  plot(lg1.mds.map.err)
  L1<-list(lg1.geno.map, 
           lg1.geno.map.err, 
           lg1.mds.map, 
           lg1.mds.map.err)
  
  names(L1) <- paste(c("geno", "geno.err", "mds", "mds.err"), i, sep = ".")
  plot_map_list(L1)
  plot_genome_vs_map(L1)
  L[[i]] <- L1
}
plot_map_list(lapply(L, function(x) x[[1]]), col = "ggstyle")


plot(L[[11]][[1]], mrk.names = T, left.lim = 40, right.lim = 70)
a1<-drop_marker(L[[11]][[1]], mrk = c("solcap_snp_c2_23987", "solcap_snp_c1_7666"))
a1<-reest_rf(a1, tol = 10e-4)
plot(a1)
a.err <- est_full_hmm_with_global_error(a1, error = 0.05)
plot(a.err)
s1 <- make_seq_mappoly(gr, arg = a[11], genomic.info = 1)
s1
tpt1 <- make_pairs_mappoly(tpt, s1)
b<-split_and_rephase(a.err, twopt = tpt1, gap.threshold = 5)
plot(b)
plot(a1)
plot(a.err)
s1 <- make_seq_mappoly(gr, arg = a[i], genomic.info = 1)
tpt1 <- make_pairs_mappoly(tpt, s1)
m1 <- make_mat_mappoly(m, s1)
##Genomic order
g.o1 <- get_genomic_order(s1)
s1.g <- make_seq_mappoly(g.o1)
plot(g.o1)
plot(m1, ord = s1.g$seq.mrk.names, fact = 3)
##3.8 minutes
system.time(lg1.geno.map<-est_rf_hmm_sequential(input.seq = s1.g,
                                                start.set = 3,
                                                thres.twopt = 10,
                                                thres.hmm = 50,
                                                extend.tail = 100, #30 for all, except 11
                                                twopt = tpt1,
                                                verbose = TRUE,
                                                tol = 10e-2,
                                                tol.final = 10e-4,
                                                phase.number.limit = 20,
                                                sub.map.size.diff.limit = 5,
                                                info.tail = TRUE,
                                                reestimate.single.ph.configuration = TRUE))
plot(lg1.geno.map)
plot(a1)
system.time(lg1.mds.map.err <- est_full_hmm_with_global_error(lg1.geno.map, error = 0.05, tol = 10e-4))
plot(lg1.mds.map.err)
plot(a.err)
plot_genome_vs_map(a.err)
plot_genome_vs_map(L[[11]][4])
