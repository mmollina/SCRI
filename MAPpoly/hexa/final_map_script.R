## Load MAPpoly
require(mappoly)

## Load screened maps
load("~/repos/SCRI/MAPpoly/hexa/final_maps.rda")

## plot maps 
# including phases
plot(final.map$ch3)
plot(final.map$ch9)
plot(final.map$ch12)

# standard view
plot_map_list(final.map, col = "ggstyle")

## Map Summary
summary_maps(final.map)

## Map distance vs. genome position
plot_genome_vs_map(final.map, same.ch.lg = TRUE)

## Compute genotype probabilities
genoprobs <- vector("list", 3)
for(i in 1:3)
  genoprobs[[i]] <- calc_genoprob_error(input.map = final.map[[i]], 
                                        step = 1, error = 0.05)
save(genoprobs, file = "~/repos/SCRI/MAPpoly/hexa/genoprob.rda")

## Homologue probabilities
h <- calc_homoprob(genoprobs)
plot(h)

## preferential pairing profiles
pp<-calc_prefpair_profiles(genoprobs)
plot(pp, type = "hom.pairs", 
     min.y.prof = 0.1, max.y.prof = 0.3, 
     P = "Beuregard", Q = "Tanzania")
