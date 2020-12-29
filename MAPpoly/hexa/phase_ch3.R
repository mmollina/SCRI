## Loading MAPpoly
require(mappoly)
load(file = "~/repos/SCRI/MAPpoly/hexa/ch3.rda")
## Phasing and re-estimating genetic map

## For hexaploids, we usually start with the following values 
# start.set = 3
# thres.twopt = 10
# thres.hmm = 10
# extend.tail = 100
# sub.map.size.diff.limit =  3

## After running the phasing for approximately 40 minutes, 
## we notice that the algorithm was testing eight linkage 
## phases using the HMM algorithm. Thus we reduced 'thres.hmm' 
## to five. We observed the same pattern. Then, we reduced it 
## to two. After that, we evaluated the gaps present in the 
## final map. This process can be very tedious, but it improves 
## the quality of the final map. You can also try to use more 
## stringent filtering in the previous steps.

system.time(lg3.map<-est_rf_hmm_sequential(input.seq = so1,
                                           start.set = 3,
                                           thres.twopt = 10,
                                           thres.hmm = 2,
                                           extend.tail = 50,
                                           twopt = tpt1,
                                           verbose = TRUE,
                                           tol = 10e-2,
                                           tol.final = 10e-3,
                                           phase.number.limit = 20,
                                           sub.map.size.diff.limit =  3,
                                           info.tail = TRUE,
                                           reestimate.single.ph.configuration = TRUE, 
                                           high.prec = TRUE))
save.image(file = "~/repos/SCRI/MAPpoly/hexa/ch3_phased.rda")
