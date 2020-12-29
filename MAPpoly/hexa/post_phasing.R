## Loading MAPpoly
require(mappoly)

## Loading data for chromosome 3
load("~/repos/SCRI/MAPpoly/hexa/ch3_phased.rda")
## Plot phased map (notice that it is a quite long map)
plot(lg3.map, phase = FALSE)

## Using the most likely map
lg3.map <- filter_map_at_hmm_thres(lg3.map, thres.hmm = 0.01)
## Obtaining the ordinary least squared map
lg3.map.ols<-reest_rf(lg3.map, input.mat = m1, tol = 10e-5, method = "ols")
## Obtaining the weighted MDS projection onto a 1D principal component
lg3.map.mds<-reest_rf(lg3.map, input.mds = o1, tol = 10e-5, method = "wMDS_to_1D_pc")
plot(lg3.map.ols)
plot(lg3.map.mds)
plot(lg3.map.mds, left.lim = 0, right.lim = 15, mrk.names = T)
plot(lg3.map.mds, left.lim = 125, right.lim = 128.7, mrk.names = T)
id<-c(lg3.map.mds$info$mrk.names[1:7], rev(lg3.map.mds$info$mrk.names)[1])
lg3.map.mds <- drop_marker(input.map = lg3.map.mds, mrk = id)
plot(lg3.map.mds)

load("~/repos/SCRI/MAPpoly/hexa/ch9_phased.rda")
plot(lg9.map, phase = F)
lg9.map.ols<-reest_rf(lg9.map, input.mat = m2, tol = 10e-5, method = "ols")
lg9.map.mds<-reest_rf(lg9.map, input.mds = o2, tol = 10e-5, method = "wMDS_to_1D_pc")
plot(lg9.map.ols)
plot(lg9.map.mds)
plot(lg9.map.mds, left.lim = 0, right.lim = 15, mrk.names = T)
plot(lg9.map.mds, left.lim = 112, right.lim = 123, mrk.names = T)
id<-c(lg9.map.mds$info$mrk.names[1:5], rev(lg9.map.mds$info$mrk.names)[1:15])
lg9.map.mds <- drop_marker(input.map = lg9.map.mds, mrk = id)
plot(lg9.map.mds)

load("~/repos/SCRI/MAPpoly/hexa/ch12_phased.rda")
plot(lg12.map, phase = F)
lg12.map.ols<-reest_rf(lg12.map, input.mat = m3, tol = 10e-5, method = "ols")
lg12.map.mds<-reest_rf(lg12.map, input.mds = o3, tol = 10e-5, method = "wMDS_to_1D_pc")
plot(lg12.map.ols)
plot(lg12.map.mds)
plot(lg12.map.mds, left.lim = 0, right.lim = 15, mrk.names = T)
plot(lg12.map.mds, left.lim = 117, right.lim = 123, mrk.names = T)
id<-c(lg12.map.mds$info$mrk.names[1:10], rev(lg12.map.mds$info$mrk.names)[1:15])
lg12.map.mds <- drop_marker(input.map = lg12.map.mds, mrk = id)
plot(lg12.map.mds)

final.map <- list(ch3 = lg3.map.mds, 
                  ch9 = lg9.map.mds, 
                  ch12 = lg12.map.mds)
save(final.map, dat, file = "~/repos/SCRI/MAPpoly/hexa/final_maps.rda")


