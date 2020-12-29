require(mappoly)
load(file = "~/repos/SCRI/MAPpoly/hexa/ch9.rda")
## Phasing and re-estimating genetic map
system.time(lg9.map<-est_rf_hmm_sequential(input.seq = so2,
                                           start.set = 3,
                                           thres.twopt = 10,
                                           thres.hmm = 10,
                                           extend.tail = 100,
                                           twopt = tpt2,
                                           verbose = TRUE,
                                           tol = 10e-2,
                                           tol.final = 10e-3,
                                           phase.number.limit = 20,
                                           sub.map.size.diff.limit =  3,
                                           info.tail = TRUE,
                                           reestimate.single.ph.configuration = TRUE, 
                                           high.prec = TRUE))
save.image(file = "~/repos/SCRI/MAPpoly/hexa/ch9_phased.rda")
