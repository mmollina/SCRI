## Load MAPpoly
require(mappoly)

## Load screened maps
load("~/repos/SCRI/MAPpoly/hexa/final_maps.rda")

## All data analysis
load("~/repos/SCRI/MAPpoly/hexa/complete_hexa_image.rda")


## plot data
plot(dat)


## plot maps 
# including phases
plot(final.map$ch3, P = "Beauregard", Q = "Tanzania")
plot(final.map$ch9, P = "Beauregard", Q = "Tanzania")
plot(final.map$ch12, P = "Beauregard", Q = "Tanzania")

# standard view
plot_map_list(final.map, col = "ggstyle")

## Map Summary
summary_maps(final.map)

## Map distance vs. genome position
plot_genome_vs_map(final.map, same.ch.lg = TRUE)

## Compute genotype probabilities
## No error modeling
genoprobs.no.err <- vector("list", 3)
for(i in 1:3)
  genoprobs.no.err[[i]] <- calc_genoprob(input.map = final.map[[i]], step = 1)
## Error modeling
genoprobs <- vector("list", 3)
for(i in 1:3)
  genoprobs[[i]] <- calc_genoprob_error(input.map = final.map[[i]], 
                                        step = 1, error = 0.05)
save(genoprobs, file = "~/repos/SCRI/MAPpoly/hexa/genoprob.rda")

## Homologue probabilities
hne <- calc_homoprob(genoprobs.no.err)
plot(hne, ind = 4)

## Homologue probabilities
h <- calc_homoprob(genoprobs[1])
plot(h, ind = 4)

## preferential pairing profiles
pp<-calc_prefpair_profiles(genoprobs)
plot(pp, type = "hom.pairs", 
     min.y.prof = 0.1, max.y.prof = 0.3, 
     P = "Beuregard", Q = "Tanzania")

updated.map <- update_map(final.map)
summary_maps(updated.map)
