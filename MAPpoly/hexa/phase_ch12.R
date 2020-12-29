require(mappoly)
load(file = "~/repos/SCRI/MAPpoly/hexa/ch12.rda")
## Phasing and re-estimating genetic map
system.time(lg12.map<-est_rf_hmm_sequential(input.seq = so3,
                                           start.set = 3,
                                           thres.twopt = 10,
                                           thres.hmm = 2,
                                           extend.tail = 50,
                                           twopt = tpt3,
                                           verbose = TRUE,
                                           tol = 10e-2,
                                           tol.final = 10e-3,
                                           phase.number.limit = 20,
                                           sub.map.size.diff.limit =  3,
                                           info.tail = TRUE,
                                           reestimate.single.ph.configuration = TRUE, 
                                           high.prec = TRUE))
save.image(file = "~/repos/SCRI/MAPpoly/hexa/ch12_phased.rda")







