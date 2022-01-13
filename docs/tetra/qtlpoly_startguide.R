## Installing and loading packages
install.packages("mappoly")
install.packages("qtlpoly")
library(mappoly)
library(qtlpoly)

## Preparing data
## Download SCRI repository: https://github.com/mmollina/SCRI
setwd("~/Documents/github/SCRI/docs/tetra/") # Workstation
setwd("~/GIT_REPOS/SCRI/docs/tetra/") # Cluster

## Loading object with conditional probabilities
## Loading genotype probabilities from MAPpoly
load("genoprobs.rda")
str(geno.prob)
geno.prob = geno.prob[c(5,7)]
## Loading phenotypes
pheno = read.table("pheno_asreml.txt", header = TRUE, row.names = 1)
head(pheno)
pheno = pheno[,1:2] # Selecting first phenotypes

## Coverting plant yield from pounds to kg
py = which(substr(colnames(pheno), start = 1, stop = 2) == "PY")
pheno[,c(py)] = pheno[,c(py)]/2.2046 

## Making sure genotype names match from pheno and geno objects
for(c in 1:length(geno.prob)) dimnames(geno.prob[[c]]$probs)[[3]] = gsub(x = dimnames(geno.prob[[c]]$probs)[[3]], pattern = ".", replacement = "-", fixed = TRUE)
all(dimnames(geno.prob[[1]]$probs)[[3]] %in% rownames(pheno))

## Importing genotypic and phenotypic data
data = read_data(ploidy = 4, geno.prob = geno.prob, pheno = pheno, step = 1)
print(data, detailed = TRUE)

## Perform random-effect multiple interval mapping (REMIM)
### Score-based resampling method to assess genome-wide significance
data.sim = simulate_qtl(data = data, mu = 0, h2.qtl = NULL, n.sim = 2, seed = 123) 
print(data.sim, detailed = TRUE)
score.null = null_model(data = data.sim$results, pheno.col = 1, n.clusters = 1)

## Running REMIM analysis
remim.mod = remim(data = data, pheno.col = 1, w.size = 20, sig.fwd = 0.20, sig.bwd = 0.05, score.null = score.null, d.sint = 1.5, n.clusters = 1)
print(remim.mod)

## Fitting final models with detected QTL
fitted.mod = fit_model(data = data, model = remim.mod, probs = "joint", polygenes = "none")
summary(fitted.mod)

## Estimating allele effects
est.effects = qtl_effects(ploidy = 4, fitted = fitted.mod)

## Saving results for VIEWpoly
save(data, file = "QTLpoly_data.RData")
save(remim.mod, file = "QTLpoly_remim.mod.RData")
save(fitted.mod, file = "QTLpoly_fitted.mod.RData") 
save(est.effects, file = "QTLpoly_est.effects.RData")

## Perform fixed-effect interval mapping (FEIM)
### Permutation method to assess genome-wide significance
perm = permutations(data = data, pheno.col = 1, n.sim = 100, n.clusters = 1)
sig.lod = numeric(length(perm$results))
for(p in 1:length(perm$results)) sig.lod[p] = quantile(sort(perm$results[[p]]), 0.95)
### Running FEIM analysis
feim.mod = feim(data = data, pheno.col = 1, w.size = 20, sig.lod = sig.lod)
print(feim.mod)
