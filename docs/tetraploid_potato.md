---
title: "Building a genetic linkage map of an autotetraploid population using MAPpoly"
author: "Marcelo Mollinari, Gabriel Gesteira, Guilherme Pereira, A Augusto Garcia, Zhao-Bang Zeng"
date: '2020-12-03'
output:
  prettydoc::html_pretty:
    theme: leonids
    highlight: github
    keep_md: yes
    toc: yes
    toc_depth: '4'
  md_document:
    variant: markdown_github
  pdf_document:
    highlight: pygments
    toc: yes
    toc_depth: '4'
  word_document:
    toc: yes
    toc_depth: '4'
linestretch: 1.2
bibliography: biblio.bib
---



# Introduction

Linkage maps are essential tools to several genetic endeavours such as quantitative trait loci (QTL) analysis, evolutionary studies, assessment of colinearity in chromosome-wise scale between genomes, and study of meiotic processes. The principle behind linkage map construction is detecting recombinant events between genomic positions and summarizing them into pairwise recombination fraction estimates. In diploids, the assessment of such a phenomenon is relatively straightforward. After the homolog duplication, sister chromatids pair up exchanging segments. The presence of informative markers, e. g. single nucleotide polymorphisms (SNPs) enable straightforward computation of the recombination fraction between pairs of genomic positions by comparing the chromosome constitution of parents and offspring. If the order of markers is unknown, it can be obtained using the pairwise recombination fractions in conjunction with optimization algorithms.

In polyploids (species with more than two sets of homolog chromosomes, or homologs), the construction of such maps is quite challenging. While in diploids a biallelic marker can be fully informative, in polyploids they allow accounting for proportions of biallelic dosages. In order to recover the multiallelic information present in polyploid species, we need to account for recombination frequencies, estimate phase configurations and reconstruct the haplotypes of both parents and individuals in the population. As the ploidy level increases, the joint computation of genomic regions becomes intensive, and other approaches such as dimension reduction need to be applied.


<!-- Building genetic linkage maps is one of the leading steps for genetic studies of any species of interest. Due to its high informativity, genetic maps provide knowledge about the physical distance between markers (independent of non-physical factors) and their linkage phases, and also allows studies about genotype-phenotype association (such as QTL mapping), study of the genetic architecture of important traits, study of evolutionary processes, helps to assemble reference genomes, and so on.

Methods to build genetic linkage maps are widespread for diploid and autotetraploid organisms. However, there is a lack of available methods for organisms with higher ploidy levels, such as sweetpotato (*Ipomoea batatas*, 6x), sugarcane (*Saccharum spp.*, 6-14x), many forage crops and other species.

Multipoint procedures, such as the Hidden Markov Chain Model (HMM), are great alternatives to estimate both linkage phases and recombination fractions between markers due to its high statistical power. Therefore, it explores information of the entire population and accounts for all possible genotype probabilities given the observed data for a linkage group. However, the dimension of these probabilities grows exponentially as ploidy increases, and for high ploidy species, procedures that avoid the need of big computational structures, such as dimensional reduction and two-point approaches, are necessary.

Building genetic linkage maps for polyploid organisms also involves estimation of their haplotypes, including both parents and progeny individuals. The haplotypes depend on recombination fractions and linkage phases, which can be estimated by both two-point and HMM approaches depending on the ploidy level. Besides genetic linkage studies, haplotypes can be very useful for QTL mapping, genomic prediction and genome wide association studies, as they allow recovering the multiallelic nature of polyploid genotypes, reduce the dimension of datasets and consequently reduce computational efforts needed to run the statistical models.

MAPpoly is an under development R package that allows building genetic linkage maps for autopolyploids with even ploidy levels. In its current version, MAPpoly can handle ploidy levels up to 8x when using the HMM approach, and up to 12x when using the two-point approach. All two-point-based functions can be run on standard computers, but we strongly recommend the use of high performance computers for HMM-based functions.-->

`MAPpoly` is an R package fully capable of building genetic linkage maps for biparental populations in polyploid species with ploidy level up to 8x. 
This package is part of the [Genomic Tools for Sweetpotato Improvement](https://sweetpotatogenomics.cals.ncsu.edu/) (GT4SP), funded by [Bill and Melinda Gates Foundation](https://www.gatesfoundation.org/). All of the procedures presented in this document are detailed in @Mollinari2019 and @Mollinari 2020. The results obtained with `MAPpoly` can be readily used for QTL mapping with the [QTLpoly package](https://github.com/guilherme-pereira/QTLpoly), which implements the procedures proposed by @Pereira2020.

Some advantages of `MAPpoly` are:

- Can handle multiple dataset types
- Can handle thousands of markers
- Does not depend on single-dose markers (SDM) to build the map
- Incorporates genomic information
- Explores multipoint information (through Hidden Markov Chain Model)
- Can handle high ploidies: up to 8x when using HMM, and up to 12x when using the two-point approach
- Can reconstruct haplotypes for parents and all individuals in the population
- Recovers the multiallelic nature of polyploid genomes
- Detects occurrence, location and frequency of multivalent pairing during meiosis
- Robust enough to build genetic linkage maps with multiallelic markers (when available)

## `MAPpoly` installation

## From CRAN (stable version)

To install MAPpoly from the The Comprehensive R Archive Network (CRAN) use

```R
install.packages("mappoly")
```

## From GitHub (development version)

You can install the development version from Git Hub. Within R, you need to install `devtools`:

```R
install.packages("devtools")
```

If you are using Windows, you must install the the latest recommended version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

To install MAPpoly from Git Hub use

```R
devtools::install_github("mmollina/mappoly", dependencies=TRUE)
```

# Loading datasets

In its current version, MAPpoly can handle the following types of datasets:

1. CSV files 
2. MAPpoly files
    - Dosage based
    - Probability based
3. [fitPoly](https://CRAN.R-project.org/package=fitPoly) files
4. VCF files

MAPpoly also is capable of importing objects generated by the following R packages 

1. [updog](https://CRAN.R-project.org/package=updog)
2. [polyRAD](https://CRAN.R-project.org/package=polyRAD)
3. [polymapR](https://CRAN.R-project.org/package=polymapR)
    - Datasets
    - Maps

Both CSV and `MAPpoly` datasets are sensible formatting errors, such as additional spaces, commas and wrong encoding (non-UTF-8). If you have any trouble, please double check your files before submitting an issue. Detailed steps of all supported files are given on the sections below. Also, a considerable proportion of redundant markers is expected on large datasets. Since markers having the same information do not provide any advantage during the mapping process, redundant ones may be removed from the dataset in order to reduce computational effort. All reading functions share the `elim.redundant` argument that automatically identifies and removes redundant markers during the analysis. The redundant markers can always be recovered and exported when the final map is ready. 

## Reading CSV files

The preparation of a CSV file for MAPpoly is relatively straightforward. It can be done in Microsoft Excel or any other spreadsheet software of your preference.  In this file, each line comprehends a marker and each column comprehends information about the marker. In its current version, `MAPpoly` can handle .csv files with allelic dosage data. 

The first line of the CSV file should contain headers for all columns. The first five columns should include the following information: marker name, dosage for both parents (one column for each), a sequence number (e.g. a chromosome number, if available) and a sequence position (e.g. the marker position within the chromosome, if available). In addition to these five headers, you should include the name of all individuals in the population. From the second line onwards, all columns should contain its values, including allelic dosages for all individuals. Missing or absent values should be represented by NA.

NOTE: If genomic information is not available, the 'sequence' and 'sequence position' columns should be filled with NA's.

Example:



**Important note: avoid spaces in .csv files.** As mentioned above, please double check your datasets for extra spaces, commas, dots and encoding. Your CSV file should be encoded using UTF-8.

You can read CSV files with the `read_geno_csv` function:


```r
ft="https://raw.githubusercontent.com/mmollina/SCRI/main/data/B2721_dose.csv"
tempfl <- tempfile()
download.file(ft, destfile = tempfl)
dat.dose.csv <- read_geno_csv(file.in  = tempfl, ploidy = 4)
```

``` r
dat.dose.csv
```

In addition to the CSV file path, you should indicate the ploidy level using the `ploidy` argument. This function automatically excludes uninformative markers. It also performs chi-square tests for all markers, considering the expected segregation patterns under Mendelian inheritance, random chromosome pairing and no double reduction. You can optionally use the `filter.non.conforming` logical argument (default = TRUE), which excludes non-expected genotypes under these assumptions. However, keep in mind that the functions in MAPpoly do not support double-reduction in its current version, and this can cause the software to abort.

## Reading MAPpoly files

Besides CSV and VCF files, MAPpoly can also handle two dataset types that follow the same format: (1) a genotype-based file (with allelic dosages) and (2) probability-based file. Both are text files with the same header, but with different genotype table formats.

For both files, the header should contain: ploidy level, number of individuals (nind), number of markers (nmrk), marker names (mrknames), individual names (indnames), allele dosages for parent 1 (dosageP), allele dosages for parent 2 (dosageQ), sequence/chromosome information (seq), position of each marker (seqpos), number of phenotypic traits (nphen) and the phenotypic data (pheno) if available. The header should be organized according to this example:

```
ploidy 4
nind 3
nmrk 5
mrknames M1 M2 M3 M4 M5
indnames Ind1 Ind2 Ind3
dosageP 0 2 0 0 3
dosageQ 1 2 1 1 3
seq 1 1 2 2 3
seqpos 100 200 50 150 80
nphen 0
pheno-----------------------
geno------------------------
```

For more information about MAPpoly file format, please see `?read_geno` and `?read_geno_prob` documentation from `MAPpoly` package.

### Using `read_geno`

The header should be followed by a table containing the genotypes (allele dosages) for each marker (rows) and for each individual (columns), as follows:

|          | Individual 1 | Individual 2 | Individual 3 |
|----------|:------------:|:------------:|:------------:|
| Marker 1 | 1            | 0            | 0            |
| Marker 2 | 3            | 0            | 2            |
| Marker 3 | 1            | 0            | 0            |
| Marker 4 | 1            | 0            | 0            |
| Marker 5 | 3            | 4            | 4            |

The final file should look like the example below:

```
ploidy 4
nind 3
nmrk 5
mrknames M1 M2 M3 M4 M5
indnames Ind1 Ind2 Ind3
dosageP 0 2 0 0 3
dosageQ 1 2 1 1 3
seq 1 1 2 2 3
seqpos 100 200 50 150 80
nphen 0
pheno-----------------------
geno------------------------
1 0 0
3 0 2
1 0 0
1 0 0
3 4 4
```

Then, use the `read_geno` function to read your file:


```r
fl = "https://raw.githubusercontent.com/mmollina/MAPpoly_vignettes/master/data/SolCAP_dosage"
tempfl <- tempfile()
download.file(fl, destfile = tempfl)
dat.dose.mpl <- read_geno(file.in  = tempfl, elim.redundant = TRUE)
```

``` r
dat.dose.mpl
```

### Using `read_geno_prob`

Following the same header described before, `read_geno_prob` reads a table containing the probability distribution for each combination of marker $\times$ individual. Each line on this table represents the combination of one marker with one individual, and its respective probabilities of having each possible allele dosage. The first two columns represent the marker and the individual, respectively, and the remaining elements represent the probability associated with each one of the possible dosages, as follows:

| Marker | Individual | $p(d=0)$ | $p(d=1)$ | $p(d=2)$ | $p(d=3)$ | $p(d=4)$ |
|--------|:----------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| M1     | Ind1       | 0.5      | 0.5      | 0.0      | 0.0      | 0.0      |
| M2     | Ind1       | 0.0      | 1.0      | 0.0      | 0.0      | 0.0      |
| M3     | Ind1       | 0.3      | 0.7      | 0.0      | 0.0      | 0.0      |
| M4     | Ind1       | 0.5      | 0.5      | 0.0      | 0.0      | 0.0      |
| M5     | Ind1       | 0.0      | 0.0      | 0.0      | 0.9      | 0.1      |
| M1     | Ind2       | 1.0      | 0.0      | 0.0      | 0.0      | 0.0      |
| M2     | Ind2       | 0.2      | 0.5      | 0.3      | 0.0      | 0.0      |
| M3     | Ind2       | 0.9      | 0.1      | 0.0      | 0.0      | 0.0      |
| M4     | Ind2       | 0.9      | 0.1      | 0.0      | 0.0      | 0.0      |
| M5     | Ind2       | 0.0      | 0.0      | 0.0      | 0.2      | 0.8      |
| M1     | Ind3       | 0.2      | 0.8      | 0.0      | 0.0      | 0.0      |
| M2     | Ind3       | 0.4      | 0.6      | 0.0      | 0.0      | 0.0      |
| M3     | Ind3       | 1.0      | 0.0      | 0.0      | 0.0      | 0.0      |
| M4     | Ind3       | 0.0      | 0.1      | 0.9      | 0.0      | 0.0      |
| M5     | Ind3       | 0.1      | 0.9      | 0.0      | 0.0      | 0.0      |

Notice that each marker $\times$ individual combination have $m+1$ associated probabilities, being $m$ the ploidy level and $m+1$ the number of possible allele dosages. Also, each line must sum to 1. The final file (header + table) should look like the following example:

```
ploidy 4
nind 3
nmrk 5
mrknames M1 M2 M3 M4 M5
indnames Ind1 Ind2 Ind3
dosageP 0 2 0 0 3
dosageQ 1 2 1 1 3
seq 1 1 2 2 3
seqpos 100 200 50 150 80
nphen 0
pheno-----------------------
geno------------------------
M1 Ind1 0.5 0.5 0.0 0.0 0.0
M2 Ind1 0.0 1.0 0.0 0.0 0.0
M3 Ind1 0.3 0.7 0.0 0.0 0.0
M4 Ind1 0.5 0.5 0.0 0.0 0.0
M5 Ind1 0.0 0.0 0.0 0.9 0.1
M1 Ind2 1.0 0.0 0.0 0.0 0.0
M2 Ind2 0.2 0.5 0.3 0.0 0.0
M3 Ind2 0.9 0.1 0.0 0.0 0.0
M4 Ind2 0.9 0.1 0.0 0.0 0.0
M5 Ind2 0.0 0.0 0.0 0.2 0.8
M1 Ind3 0.2 0.8 0.0 0.0 0.0
M2 Ind3 0.4 0.6 0.0 0.0 0.0
M3 Ind3 1.0 0.0 0.0 0.0 0.0
M4 Ind3 0.0 0.1 0.9 0.0 0.0
M5 Ind3 0.1 0.9 0.0 0.0 0.0
```

To read the dataset, one should use:


```r
ft="https://raw.githubusercontent.com/mmollina/MAPpoly_vignettes/master/data/SolCAP"
tempfl <- tempfile()
download.file(ft, destfile = tempfl)
dat.prob.mpl <- read_geno_prob(file.in  = tempfl, prob.thres = 0.95, elim.redundant = TRUE)
```

``` r
dat.prob.mpl
```

**Important note:** as this type of file contains the probability distribution of all possible dosages, it will take longer to read. 

This function automatically excludes monomorphic markers, keeping only informative ones. You can define the minimum probability value necessary to call a dosage using the `prob.thres` argument. If the higher probability for a marker $\times$ individual passes this threshold, then its associated dosage is used. However, if none of the probabilities reach this threshold, then its dosage is considered missing (NA).
This function also performs chi-square tests for all markers, considering the expected segregation patterns under Mendelian inheritance, random chromosome pairing and no double reduction. You can optionally use the `filter.non.conforming` logical argument (default = TRUE), which excludes non-expected genotypes under these assumptions.

## Reading VCF files {#read_vcf}

VCF files are less sensible to errors, because they are usually produced by automated SNP calling pipelines and less susceptible to user edition. `MAPpoly` can also handle VCF files of version 4.0 and higher produced by the most common softwares, such as TASSEL, GATK, Stacks and many others. As few of these softwares can handle poliploidy and estimate genotypes correctly, you may use other softwares that are dedicated to estimate the allele dosages. Briefly, these softwares model the ratio between allele read counts (or intensity) for each marker $\times$ individual combination, and determines which is the most probable allele dosage given the observed ratio and other *a priori* information. Examples of these softwares are [SuperMASSA](http://statgen.esalq.usp.br/SuperMASSA/), [fitPoly](https://cran.r-project.org/web/packages/fitPoly/index.html), [ClusterCall](https://potatobreeding.cals.wisc.edu/wp-content/uploads/sites/161/2017/08/ClusterCall_Download.zip), [updog](https://cloud.r-project.org/web/packages/updog/index.html), [PolyRAD](https://cran.r-project.org/web/packages/polyRAD/vignettes/polyRADtutorial.html) and many others. After allele dosage estimation, your VCF file should contain GT values like **1/1/1/0** (for an autotetraploid) rather than **1/0**. Since `MAPpoly` uses allelic dosages (or their probabilities) to build genetic maps, we strongly recommend that you use one of these softwares to estimate allele dosages before building the map. Both `PolyRAD` and `updog` have direct integration with `MAPpoly`, as described in the next section.

To demonstrate the `read_vcf` function, lets download an autohexaploid sweetpotato VCF file from the MAPpoly's vignettes repository on Github and read it:


```r
download.file("https://github.com/mmollina/MAPpoly_vignettes/raw/master/data/BT/sweetpotato_chr1.vcf.gz", destfile = 'chr1.vcf.gz')
dat.dose.vcf = read_vcf(file = 'chr1.vcf.gz', parent.1 = "PARENT1", parent.2 = "PARENT2")
```

``` r
dat.dose.vcf
```

Besides the path to your VCF file, you should indicate `parent.1` and `parent.2` names. Please notice that their names must be exactly the same strings that appear in your VCF file. The ploidy level will be automatically detected, but you may indicate it using the optional `ploidy` argument to let the function check for possible errors. For species with variable ploidy levels (i.e. sugarcane), please indicate the desired ploidy level using the `ploidy` argument; if absent, `MAPpoly` will use the smallest ploidy level detected. 

This function also has options to filter out undesired markers or data points, such as those with low depth or high proportion of missing data. You can define the following filter arguments: set `min.av.depth` to any integer level in order to remove markers that show average depth below this value (default = 0); set `min.gt.depth` to any integer level in order to remove data points that present depth below this value (default = 0); set `max.missing` to any value between 0 and 1 in order to remove markers that present missing data proportion above this value (default = 1).

`read_vcf` performs chi-square tests for all markers, considering the expected segregation patterns under Mendelian inheritance, random chromosome pairing and no double reduction. You can optionally use the `filter.non.conforming` logical argument (default = TRUE), which excludes non-expected genotypes under these assumptions. The p-value threshold used by the segregation test can be defined by the `thresh.line` argument.

Please notice that the returning object from `read_vcf` has some additional information: reference and alternative alleles (bases) for each marker; and average depth of each marker. You can inspect all marker depths using the following code as an example:


```r
library(ggplot2)
dosage_combs = cbind(dat.dose.vcf$dosage.p, dat.dose.vcf$dosage.q)
dc_simplex = apply(dosage_combs,1,function(x) if(all(c(0,1) %in% x) | all(c(dat.dose.vcf$m-1, dat.dose.vcf$m) %in% x)) return(TRUE) else return(FALSE))
dc_dsimplex = apply(dosage_combs,1,function(x) if(all(x == c(1,1)) | all(x == c(dat.dose.vcf$m-1, dat.dose.vcf$m-1))) return(TRUE) else return(FALSE))

dc_simplex[which(dc_simplex == TRUE)] = "simplex"
dc_simplex[which(dc_dsimplex == TRUE)] = 'double simplex'
dc_simplex[which(dc_simplex == 'FALSE')] = 'multiplex'

data_depths = data.frame('Marker depths' = dat.dose.vcf$all.mrk.depth,
                         'Depth classes' = findInterval(dat.dose.vcf$all.mrk.depth, c(200,300,400,500,600,50000)),
                         'Dosage combinations' = dc_simplex, check.names = F)

ggplot(data_depths, aes(fill=`Dosage combinations`, x=`Depth classes`, y=`Marker depths`)) +
  geom_bar(position = 'stack', stat = 'identity') +
  scale_x_continuous(breaks=0:5, labels=c(">200","200-300","300-400","400-500","500-600", ">600"))
```


## Importing data from third party packages

As mentioned above, polyploid sequencing data must be used to estimate allelic dosages prior to linkage map building. The R package `PolyRAD` has its own function to export genotypes to the `MAPpoly`'s genotype probability distribution format. One may use the commands above to import from `PolyRAD`:


```r
# load example dataset from polyRAD
library(polyRAD)
data(exampleRAD_mapping)
exampleRAD_mapping = SetDonorParent(exampleRAD_mapping, "parent1")
exampleRAD_mapping = SetRecurrentParent(exampleRAD_mapping, "parent2")
exampleRAD_mapping = PipelineMapping2Parents(exampleRAD_mapping)

# export to MAPpoly
outfile2 = tempfile()
Export_MAPpoly(exampleRAD_mapping, file = outfile2)

# Read in MAPpoly
mydata_polyrad = read_geno_prob(outfile2)
```

``` r
mydata_polyrad
```

You can also use the `MAPpoly`'s function `import_from_updog` to import any dataset generated by `updog`'s function `multidog`, following the example below:


```r
# Load example dataset from updog
library(updog)
data(uitdewilligen)
mout = multidog(refmat = t(uitdewilligen$refmat), 
                sizemat = t(uitdewilligen$sizemat), 
                ploidy = uitdewilligen$ploidy, 
                model = "f1",
                p1_id = colnames(t(uitdewilligen$sizemat))[1],
                p2_id = colnames(t(uitdewilligen$sizemat))[2],
                nc = 4)
mydata_updog = import_from_updog(mout)
```

``` r
mydata_updog
```

Please notice that `updog` removes both sequence and sequence position information that may be present in the VCF file. We highly recommend that you use this information during the linkage map building, when available.

## Combining multiple datasets

It is not rare to have two or more datasets regarding the same population or individuals from different sources of molecular data, such as SNP chips, GBS and/or microsatellites. `MAPpoly` can combine datasets when individual names are the same, using the function `merge_datasets`. To demonstrate its functionality, lets download two VCF files (autohexaploid sweetpotato) from the MAPpoly vignettes repository on Github, and read them using `read_vcf` function:


```r
# Downloading VCF files regarding chromosome 1 and 2
download.file("https://github.com/mmollina/MAPpoly_vignettes/raw/master/data/BT/sweetpotato_chr1.vcf.gz", destfile = 'chr1.vcf.gz')
download.file("https://github.com/mmollina/MAPpoly_vignettes/raw/master/data/BT/sweetpotato_chr2.vcf.gz", destfile = 'chr2.vcf.gz')
data1 = read_vcf(file = 'chr1.vcf.gz', parent.1 = "PARENT1", parent.2 = "PARENT2")
data2 = read_vcf(file = 'chr2.vcf.gz', parent.1 = "PARENT1", parent.2 = "PARENT2")
```

As we can see, both have different markers for the same population:


```r
# See datasets
print(data1)
print(data2)
```

Now lets merge them and see the output:


```r
# Merge datasets
merged_data = merge_datasets(data1, data2)
```


```r
print(merged_data)
```

Notice that all markers of both datasets were merged successfully, which allows using just one (merged) dataset in the following steps of map construction. 

## Exploratory Analysis 

### Whole dataset

For didatic purposes, we will keep using the tetraploid potato array data (loaded using the examples above). We will construct a genetic map of the B2721 population, which is a cross between two tetraploid potato varieties: Atlantic and B1829-5. The population comprises 160 offsprings genotyped with the SolCAP Infinium 8303 potato array. The dataset also contains the genomic order of the SNPs from the _Solanum tuberosum_ genome version 4.03. The genotype calling was performed with fitPoly R package using [this pipeline](https://github.com/mmollina/Autopolyploid_Linkage/blob/master/src/solcap_map_construction/snp_calling/genotype_calling_public_data_fittetra.R). Another option would be to use ClusterCall and [this pipeline](https://mmollina.github.io/tutorials/solcap/solcap_example.html).

Once the data is loaded, you can explore the dataset using the `print` function:


```r
print(dat.dose.mpl, detailed = TRUE)
```

This function outputs information about the dataset including the ploidy level, total number of individuals, total number of markers, number of informative markers, proportion of missing data and redundant markers. If `detailed = TRUE`, the function also outputs the number of markers in each sequence, if available, and the number of markers contained in all possible dosage combinations between both parents.

You can also explore the dataset visually using the `plot` function:


```r
plot(dat.dose.mpl)
```

The output figure shows a bar plot on the left-hand side with the number of markers contained in each allele dosage combination between both parents. The right labels indicate allele dosages for Parent 1 and Parent 2, respectively. The upper-right plot contains the $\log_{10}(p-value)$ from $\chi^2$ tests for all markers, considering the expected segregation patterns under Mendelian inheritance. The lower-right plot contains a graphical representation of the allele dosages and missing data distribution for all markers and individuals. Finally, the bottom-right graphic shows the proportion of redundant markers in the dataset, when available.

### Marker-specific

If you want to view a specific marker information, use the `plot_mrk_info` function. You should indicate your dataset object using the `input.data` argument, and the desired marker using the `mrk` argument. You can indicate the marker using its number or its name (string):


```r
# For a dosage-based data analysis of marker 104
plot_mrk_info(input.data = dat.dose.mpl, mrk = 240)
```


```r
# For a probability-based data analysis of the marker solcap_snp_c1_13686
plot_mrk_info(input.data = dat.prob.mpl, mrk = 'solcap_snp_c1_13686')
```

When applied to a dosage-based dataset, the function outputs a figure showing: marker name and position in the dataset, allele dosage in parents 1 and 2, proportion of missing data, p-value of the associated $\chi^2$ test for Mendelian segregation, sequence and position information (when available). The figure also contains a plot with the allele dosage and missing data distibution in the population.

When applied to a probability-baed dataset, the function outputs the probability threshold and a
3-dimensional plot containing the probability distribution for each allele dosage considering all individuals.

# Filtering and Quality Control

In order to build a good genetic map, good quality data is desired to guarantee reliable estimates of recombination fractions, linkage phases and haplotypes. High proportions of messy data, such as unexpected segregation patterns and missing data, will reduce the reliability of these estimates. Furthermore, sequencing technologies are able to produce hundreds of thousands of markers, and a considerable proportion of redundant markers is expected. Markers that carry the same information are not informative and must be removed from the dataset, in order to reduce computational effort. `MAPpoly` handle some filtering functions that are described in the sections below.

## Missing data filtering

`MAPpoly` is able to filter markers and/or individuals that exceeds a defined threshold for missing data. The function `filter_missing` does that and creates or updates your dataset according to its arguments' values. The argument `input.data` should contain your dataset object, and you can choose to filter either by 'marker' or 'individual' using the `type` argument (string). You can also define the maximum proportion of missing data using the `filter.thres` argument (ranges from 0 to 1, i.e. a threshold of 0.2 will keep just markers or individuals with less than 20% of missing data). When TRUE (default), the `inter` argument plots markers or individuals vs. frequency of missing data.


```r
# Filtering dataset by marker
dat.filt.mrk <- filter_missing(input.data = dat.dose.mpl, type = "marker", 
                               filter.thres = 0.2, inter = TRUE)
```


```r
print(dat.filt.mrk)
```


```r
# Filtering dataset by individual
dat.filt.ind <- filter_missing(input.data = dat.filt.mrk, type = "individual", 
                               filter.thres = 0.1, inter = TRUE)                            
```


```r
print(dat.filt.ind)
```


```r
dat.dose.filt <- dat.filt.ind
```

In this dataset, just 43 markers presented a proportion of missing data above the defined threshold, while 16 individuals exceeded the defined threshold. Then we will use the final filtered dataset during the rest of the analysis. Please notice that the function `read_vcf` also provides parameters to filter out markers depending on their average depth and missing data, as well as removes data points that do not reach a minimum pre-defined depth value. Check the [`read_vcf`](#read_vcf) section for more information.

## Segregation test

Another very important point to consider is the expected marker segregation pattern under Mendelian inheritance. Markers with messy or distorted segregation can produce unreliable estimates and may be removed (at least for a while) from the dataset. A good test for that is the chi-square ($\chi^2$), which basically matches expected genotype frequencies against observed frequencies and calculates the associated p-value. In order to define the p-value threshold for the tests, we will use the Bonferroni correction: 

$$\alpha_{thres} = \frac{\alpha}{\#markers}$$

We will also assume that only random chromosome bivalent pairing occurs and there is no double reduction.


```r
pval.bonf <- 0.05/dat.dose.filt$n.mrk
mrks.chi.filt <- filter_segregation(dat.dose.filt, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init <- make_seq_mappoly(mrks.chi.filt)
```

Please notice that `filter_segregation` does not produce a filtered dataset; it just tells you which markers follow the expected Mendelian segregation pattern. To select these markers from your dataset, you may use the `make_seq_mappoly` function.


```r
plot(seq.init)
```

It is worth to mention that all redundant markers identified during data reading step are stored in the main dataset. The redundant markers are automatically removed and can be added back once the maps are finished using the function `update_map` (described in details later during this tutorial).

# Two-point analysis

Once the markers are filtered and selected, we need to compute the pairwise recombination fraction between all of them (two-point analysis). First, let us load the genotype counts ($\zeta_{\mbox{T}_{k},\mbox{T}_{k^{\prime}}}(l_{P}, l_{Q})$) defined in equation 20 in [Mollinari and Garcia, 2019](https://doi.org/10.1534/g3.119.400378). This object is fundamental to perform the dimension reduction of the transition space.


```r
counts <- cache_counts_twopt(input.seq = seq.init, cached = TRUE)
```

``` r
counts
```

When `cached = TRUE`, the function loads all counts from a internal file instead of performing all calculations again. 

Next, the function `est_pairwise_rf` is used to estimate all the pairwise recombination fractions between markers in the sequence provided. Since the output object is too big to be fully displayed on the screen, `MAPpoly` shows a summary. Please notice that parallel computation is available and, if you want to use of that, you need to define the available number of cores in your machine and also guarantee that you have sufficient RAM memory for that. Remember that it is always good to leave one core available for the system, and be aware that this step **will take a while** to compute if you have few cores available.


```r
# Defining number of cores
n.cores = parallel::detectCores() - 1

#(~ 9.5 minutes using 14 cores)
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, 
                                   count.cache = counts, 
                                   n.clusters = n.cores)
```


```r
all.rf.pairwise
```

To assess the recombination fraction between a particular pair of markers, say markers 93 and 98, we use the following syntax:


```r
all.rf.pairwise$pairwise$`93-98`
plot(all.rf.pairwise, first.mrk = 93, second.mrk = 98)
```

In this case, `93-98` represents the position of the markers in the filtered data set. The name of the rows have the form `x-y`, where `x` and `y` indicate how many homologous chromosomes share the same allelic variant in parents $P1$ and $P2$, respectively (see [Mollinari and Garcia, 2019](https://doi.org/10.1534/g3.119.400378) for notation). The first column indicates the LOD Score in relation to the most likely linkage phase configuration. The second column shows the estimated recombination fraction for each configuration, and the third indicates the LOD Score comparing the likelihood under no linkage ($r = 0.5$) with the estimated recombination fraction (evidence of linkage).

## Assembling recombination fraction and LOD Score matrices

Recombination fraction and LOD Score matrices are fundamental in genetic mapping. Later in this tutorial, we will use these matrices as the basic information to order markers and also to perform some diagnostics. To convert the two-point object into recombination fraction and LOD Score matrices, we need to assume thresholds for the three columns observed in the previous output. The arguments `thresh.LOD.ph` and `thresh.LOD.rf` set LOD Scores thresholds for the second most likely linkage phase configuration and recombination fraction. Here we assume `thresh.LOD.ph = 0` and `thresh.LOD.rf = 0`, thus no matter how likely is the second best option, all the computed values will be considered. The argument `thresh.rf = 0.5` indicates that the maximum accepted recombination fraction is `0.5`. To convert these values in a recombination fraction matrix, we use the function `rf_list_to_matrix`.


```r
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise, n.clusters = n.cores)
```

Please notice that for big datasets, you may use the multi-core support to perform the conversions, using the parameter `n.clusters` to define the number of CPU's. It is also possible to filter again for the parameters mentioned above, such as `thresh.LOD.ph`, `thresh.LOD.rf` and `thresh.rf`.

We can also plot this matrix using the reference genome order. For doing so, we use the function `get_genomic_order` to get the genomic order of the input sequence and use the resulting order to index the recombination fraction matrix. If the reference order is consistent with the marker order in this specific population, we should observe a block-diagonal matrix and within each sub-matrix, a monotonic pattern. 


```r
id<-get_genomic_order(seq.init)
plot(mat, ord = rownames(id), index = FALSE)
```



As expected, we can observe the block-diagonal and monotonic patterns. In the previous case, the thresholds allowed to plot almost all points in the recombination fraction matrix. The empty cells in the matrix indicate markers where it is impossible to detect recombinant events using two-point estimates (e.g., between $1 \times 0$ and $0 \times 1$ marker). Yet, if the thresholds become more stringent (higher LODs and lower rf), the matrix becomes more sparse.


# Assembling linkage groups

The function `group_mappoly` assign markers to linkage groups using the recombination fraction matrix obtained above. The user can provide an expected number of groups or run the interactive version of the function using `inter = TRUE`. Since in this data set we expect 12 linkage groups (basic chromosome number in potato), we use `expected.groups = 12`. If the data set provides the chromosome where the markers are located, the function allows comparing the groups obtained using the pairwise recombination fraction and the chromosome information provided using the `comp.mat = TRUE`.


```r
grs <- group_mappoly(input.mat = mat,
                     expected.groups = 12,
                     comp.mat = TRUE, 
                     inter = TRUE)
```

``` r
grs
```

Here, we have the 3639 markers distributed in 12 linkage groups. The rows indicate linkage groups obtained using linkage information and the columns are the chromosomes in the reference genome. Notice the diagonal indicating the concordance between the two sources of information. Now, we can plot the resulting marker cluster analysis.


```r
plot(grs)
```

Once the linkage groups are properly assembled, we use the function `make_seq_mappoly` to make marker sequences from the group analysis. We will assemble a list with 12 positions, each one containing the corresponding linkage group sequence. Also, we will use only markers allocated in the diagonal of the previous comparison matrix. Thus only markers that were assigned to a particular linkage group using both sources of information will be considered. We will also assemble a smaller two-point object using the functions `make_pairs_mappoly` and `rf_snp_filter` to facilitate further parallelization procedures.


```r
LGS<-vector("list", 12)
for(j in 1:12){
    temp1 <- make_seq_mappoly(grs, j)
    temp2 <- get_genomic_order(temp1) # assembling sequence considering the genomic order
    lg.id <- as.numeric(names(which.max(table(temp2[,1]))))
    nm <- rownames(temp2)[which(temp2[,1] == lg.id)]
    temp3 <- make_seq_mappoly(dat.dose.filt, nm)
    tpt <- make_pairs_mappoly(all.rf.pairwise, input.seq = temp3)
    lgtemp <- rf_snp_filter(input.twopt = tpt)
    LGS[[lg.id]] <- list(lg = lgtemp, 
        tpt = make_pairs_mappoly(all.rf.pairwise, input.seq = lgtemp))
}
```

Now, let us print the recombination fraction matrices for each linkage group.



# Estimating the map for a given order

In this section, we will use the marker order provided by the _Solanum tuberosum_ genome version 4.03. The MDS ordering approach will be addressed later in this tutorial. The estimation of the genetic map for a given order involves the computation of recombination fraction between adjacent markers and also finding the linkage phase configuration of those markers in both parents. The core function to perform these tasks in `MAPpoly` is `est_rf_hmm_sequential`. This function uses the pairwise recombination fraction as the first source of information to sequentially position allelic variants in specific homologs. For situations where pairwise analysis has limited power, the algorithm relies on the likelihood obtained through a hidden Markov model (HMM) [Mollinari and Garcia, 2019](https://doi.org/10.1534/g3.119.400378). Once all markers are positioned, the final map is reconstructed using the HMM multipoint algorithm. 


Several arguments are available to control the inclusion and phasing of the markers in the chain. The argument `start.set` defines the number of initial markers to be used in a exhaustive search for the most probable configuration. `thres.twopt` receives the threshold to whether when the linkage phases compared via two-point analysis should be considered, and the HMM analysis should not be used to infer the linkage phase (A. K. A. $\eta$ in [Mollinari and Garcia, 2019](https://doi.org/10.1534/g3.119.400378)). `thres.hmm` receives the threshold for keeping competing maps computed using HMM (if the two-point analysis was not enough) in the next round of marker insertion. `extend.tail` indicates the number of markers that should be considered at the end of the chain to insert a new marker. `tol` and `tol.final` receive the desired accuracy to estimate the sub-maps during the sequential phasing procedure and the desired accuracy in the final map. `phase.number.limit` receives the limit number of linkage phase configurations to be tested using HMM. `info.tail` is a logical argument: if TRUE it uses the complete informative tail (last markers in the chain that allow all homologous to be distinguished in the parents) of the chain to calculate the likelihood of the linkage phases. 

First, as an example, let us estimate the map for linkage group 3. The values used for all arguments were obtained using a balance of processing speed and accuracy of the algorithm. As an exercise, it is interesting to try different values and check out the results. For now, let us stick with the following values (**this step can take some hours to finish**, depending on chosen parameters):


```r
lg3.map <- est_rf_hmm_sequential(input.seq = LGS[[3]]$lg,
                                start.set = 10,
                                thres.twopt = 10, 
                                thres.hmm = 10,
                                extend.tail = 200,
                                info.tail = TRUE, 
                                twopt = LGS[[3]]$tpt,
                                sub.map.size.diff.limit = 10, 
                                phase.number.limit = 20,
                                reestimate.single.ph.configuration = TRUE,
                                tol = 10e-3,
                                tol.final = 10e-4)
```

Now, use the functions `print` and `plot` to view the results:


```r
print(lg3.map)
plot(lg3.map)
```

Colored rectangles (red and blue) indicate the presence of the allelic variant in each one of the four homologous in both parents, $P_1$ and $P_2$.

## Reestimating the map considering genotyping errors

Though current technologies enabled the genotyping of thousands of SNPs, they are quite prone to genotyping errors, especially in polyploid species. One way to address this problem is to associate a probability distribution to each one of the markers and allow the HMM to update their probability. This procedure can be applied using either the probability distribution provided by the genotype calling software (loaded in `MAPpoly` using the function `read_geno_prob`) or assuming a global genotype error. For a detailed explanation of this procedure, please see [Mollinari and Garcia, 2019](https://doi.org/10.1534/g3.119.400378). Briefly, the use of a prior information will update the genotype of the markers based on a global chromosome structure. In this tutorial, since we are analyzing the dosage data with no probability distribution associated, the second approch will be used. 


```r
lg3.map.error <- est_full_hmm_with_global_error(input.map = lg3.map, error = 0.05)
```


```r
plot(lg3.map.error)
```

Notice that a global genotyping error of 5% was used, and the resulting map was smaller than the previous one. Also, some markers were "attracted" and some markers were "repealed" as a result of the smaller confidence used for each marker genotype.

## Reinserting redundant markers

As mentioned before, redundant markers are automatically removed during the analysis to reduce computational calculations. Besides that, one may want to see the genetic linkage maps considering all markers in the full dataset, including the redundant ones. The addition of redundant markers to a map can be done by just calling the function `update_map`. Here we will show the addition in the previous map:



```r
lg3.map.updated = update_map(lg3.map)
lg3.map
lg3.map.updated
```

As can be seen, both maps are identical except for the number of markers.

# Ordering markers using MDS and reestimating the map

So far the map was reestimated using the genomic order. In real situations, unless a genomic information is provided,  the markers need to be ordered using an optimization technique. Here, we use the MDS (multidimensional scaling) algorithm, proposed in the context of genetic mapping by @Preedy2016. The MDS algorithm requires a recombination fraction matrix, which will be transformed in distances using a mapping function (in this case the Haldane's mapping function). First, let us gather the pairwise recombination fractions for all three linkage groups:


```r
mt <- lapply(LGS, function(x) rf_list_to_matrix(x$tpt))
```

Now, for each matrix contained in the object `mt`, we use the MDS algorithm:


```r
mds.ord <- lapply(mt, mds_mappoly)
```

Usually at this point, the user can make use of diagnostic plots to remove markers that are disturbing the ordering procedure. Here we didn't use that procedure, but we encourage the user to check the example in `?mds_mappoly`. Now, let us compare the estimated and the genomic orders (feel free to run the last commented line and see interactive plots):


```r
LGS.mds<-vector("list", 12)
for(j in 1:12){
  lgtemp <- make_seq_mappoly(mds.ord[[j]])
  LGS.mds[[j]] <- list(lg = lgtemp, 
      tpt = make_pairs_mappoly(all.rf.pairwise, input.seq = lgtemp))
}
```


```r
geno.vs.mds <- NULL
for(i in 1:length(LGS.mds)){
  geno.vs.mds<-rbind(geno.vs.mds,
                     data.frame(mrk.names = LGS.mds[[i]]$lg$seq.mrk.names,
                                mds.pos = seq_along(LGS.mds[[i]]$lg$seq.mrk.names),
                                genomic.pos = order(LGS.mds[[i]]$lg$sequence.pos),
                                LG = paste0("LG_", i)))
}

require(ggplot2)
p<-ggplot(geno.vs.mds, aes(genomic.pos, mds.pos)) +
  geom_point(alpha = 1/5, aes(colour = LG)) +
  facet_wrap(~LG) +  xlab("Genome Order") + ylab("Map Order")
p
#plotly::ggplotly(p)
```

Although several local inconsistencies occurred, the global diagonal pattern indicates a consistent order for all linkage groups using both approaches. Now, let us build the genetic map of linkage group 3 using the MDS order (remember that **this can take a while to finish**):


```r
lg3.map.mds <- est_rf_hmm_sequential(input.seq = LGS.mds[[3]]$lg,
                                start.set = 10,
                                thres.twopt = 10, 
                                thres.hmm = 10,
                                extend.tail = 200,
                                info.tail = TRUE, 
                                twopt = LGS.mds[[3]]$tpt,
                                sub.map.size.diff.limit = 10, 
                                phase.number.limit = 20,
                                reestimate.single.ph.configuration = TRUE,
                                tol = 10e-3,
                                tol.final = 10e-4)
```

And plot the map:


```r
plot(lg3.map.mds)
```

It is also possible to compare the maps using both genomic-based and MDS-based orders with the function `plot_map_list`:


```r
plot_map_list(list(lg3.map, lg3.map.mds), col = c("red", "blue"), title = "")
```

The genomic-based map included the same number of markers but is smaller than the MDS-based map, which indicates a better result. To formally compare the maps, one needs to select the markers that are present in both maps. Interestingly enough, both procedures included the same markers in the final map; however, we provide the code to perform the comparison even if the maps share only a subset of markers:


```r
mrks.in.gen<-intersect(lg3.map$maps[[1]]$seq.num, lg3.map.mds$maps[[1]]$seq.num)
mrks.in.mds<-intersect(lg3.map.mds$maps[[1]]$seq.num, lg3.map$maps[[1]]$seq.num)
if(cor(mrks.in.gen, mrks.in.mds) < 0){
  mrks.in.mds <- rev(mrks.in.mds)
  lg3.map.mds <- rev_map(lg3.map.mds)
}
map.comp.3.gen<-get_submap(input.map = lg3.map, match(mrks.in.gen, lg3.map$maps[[1]]$seq.num), verbose = FALSE)
map.comp.3.mds<-get_submap(input.map = lg3.map.mds, match(mrks.in.mds, lg3.map.mds$maps[[1]]$seq.num), verbose = FALSE)
prob.3.gen<-extract_map(lg3.map)
prob.3.mds<-extract_map(lg3.map.mds)
names(prob.3.gen)<-map.comp.3.gen$maps[[1]]$seq.num
names(prob.3.mds)<-map.comp.3.mds$maps[[1]]$seq.num
```


```r
matplot(t(data.frame(prob.3.gen,prob.3.mds[names(prob.3.gen)])), 
        type="b", pch="_", col=1, lty=1, lwd = .5, xlab= "", 
        ylab="Marker position (cM)", axes = F)
axis(2)
mtext(text = round(map.comp.3.gen$maps[[1]]$loglike,1), side = 1, adj = 0)
mtext(text = round(map.comp.3.mds$maps[[1]]$loglike,1), side = 1, adj = 1)
mtext(text = "Genomic", side = 3, adj = 0)
mtext(text = "MDS", side = 3, adj = 1)
```

Please notice that these maps have the same local inversions shown in the dot plots presented earlier. In this case, the log-likelihood of the genomic order is higher than the one obtained using the MDS order. For this linkage group, the genomic-based map was chosen as the best one.

# Parallel map construction

## Using one core by LG

Now, the mapping procedure will be applied to all linkage groups using parallelization. Although users must be encouraged to compare both MDS and genomic orders following the previous example, the genomic order will be considered as an example (remember that this step will take a **long time** to run):


```r
## Performing parallel computation
my.phase.func<-function(X){
  x<-est_rf_hmm_sequential(input.seq = X$lg,
                                start.set = 10,
                                thres.twopt = 10, 
                                thres.hmm = 10,
                                extend.tail = 200,
                                info.tail = TRUE, 
                                twopt = X$tpt,
                                sub.map.size.diff.limit = 8, 
                                phase.number.limit = 10,
                                reestimate.single.ph.configuration = TRUE,
                                tol = 10e-3,
                                tol.final = 10e-4)
  return(x)
}
system.time({
  cl <- parallel::makeCluster(n.cores)
  parallel::clusterEvalQ(cl, require(mappoly))
  parallel::clusterExport(cl, "dat.dose.filt")
  MAPs <- parallel::parLapply(cl,LGS,my.phase.func)
  parallel::stopCluster(cl)
})
```

A traditional linkage map plot can be generated including all linkage groups, using the function `plot_map_list`:


```r
plot_map_list(MAPs, col = "ggstyle")
```

Following the reconstruction of LG 3 shown before, let us consider a global genotyping error of 5% to reestimate the final maps:


```r
my.error.func<-function(X){
  x<-est_full_hmm_with_global_error(input.map = X, 
                                    error = 0.05, 
                                    tol = 10e-4, 
                                    verbose = FALSE)
  return(x)
}
system.time({
  cl <- parallel::makeCluster(n.cores)
  parallel::clusterEvalQ(cl, require(mappoly))
  parallel::clusterExport(cl, "dat.dose.filt")
  MAP.err <- parallel::parLapply(cl,MAPs,my.error.func)
  parallel::stopCluster(cl)
})
```

Comparing both results:


```r
all.MAPs <- NULL
for(i in 1:12) 
  all.MAPs<-c(all.MAPs, MAPs[i], MAP.err[i])
plot_map_list(map.list = all.MAPs, col = rep(c("#E69F00", "#56B4E9"), 12))
```

Then, the map that included modeling of genotype errors was chosen as the best one.


```r
plot_map_list(MAP.err)
```

<!-- ## Using submap features -->

<!-- In the previous section, the parallelization procedure of one chromosome by core was demonstrated. Now, another way to build the genetic linkage map will be shown: based on a known sequence, the map is divided in batches (or submaps) and all steps are performed for these batches. After all calculations, the submaps are joined again to reconstruct the full map. This can be very useful to speed up the map building process based on a known ordered sequence. Here, the method applied to the chromosome 3 of the same tetraploid potato dataset will be shown. You may adjust the number of cores if you want to try this using your personal machine. -->

<!-- ```{r batches_definition} -->
<!-- ## Subsetting chromosome 3 -->
<!-- seq_submaps = make_seq_mappoly(dat.dose.filt, 'seq3') -->
<!-- markers_submaps = intersect(seq_submaps$seq.num, all.rf.pairwise$seq.num) -->
<!-- seq_submaps = make_seq_mappoly(dat.dose.filt, markers_submaps) -->
<!-- tpt_submaps = make_pairs_mappoly(all.rf.pairwise, input.seq = seq_submaps) -->
<!-- seq_submaps_filter = rf_snp_filter(input.twopt = tpt_submaps) -->
<!-- markers_submaps = seq_submaps_filter$seq.mrk.names -->

<!-- ## Running in parallel with submaps -->
<!-- submaps_merged = est_map_parallel(data = dat.dose.filt,  -->
<!--                                   markers = markers_submaps,  -->
<!--                                   partial_tpt = tpt_submaps,  -->
<!--                                   n.batches = 10, -->
<!--                                   n.cores = 10, -->
<!--                                   start.set = 5, -->
<!--                                   thres.twopt = 1,  -->
<!--                                   thres.hmm = 1,  -->
<!--                                   submap.size.diff = 3, -->
<!--                                   thres.twopt2 = 2,  -->
<!--                                   thres.hmm2 = 2, -->
<!--                                   tol = 1e-1) -->
<!-- ``` -->

<!-- Now lets compare the number of markers in the final maps: -->

<!-- ```{r, eval=T, fig.width=10} -->
<!-- lg3.map -->
<!-- submaps_merged[[1]] -->
<!-- same_markers = as.character(intersect(submaps_merged[[1]]$maps[[1]]$seq.num, lg3.map$maps[[1]]$seq.num)) -->
<!-- length(same_markers) -->
<!-- ``` -->
<!-- The first map (built with sequential approach) contained 259 markers, while the second (built using submaps feature) contained 255 markers, which of 247 are common to both maps. Considering that both were constructed using different approaches, this difference is expected. Now lets see the second map: -->

<!-- ```{r, eval=T, fig.width=10} -->
<!-- plot(submaps_merged[[1]]) -->
<!-- ``` -->

<!-- We can visually compare both maps: -->

<!-- ```{r, eval=T, fig.width=10} -->
<!-- maps = list(geno.map = lg3.map, geno.submaps = submaps_merged[[1]]) -->
<!-- plot_map_list(maps, col = "ggstyle") -->
<!-- ``` -->

<!-- And check for phase configurations between both maps: -->

<!-- ```{r, eval=T, fig.width=10} -->
<!-- plot_compare_haplotypes(4, -->
<!--                         lg3.map$maps[[1]]$seq.ph$P[same_markers], -->
<!--                         lg3.map$maps[[1]]$seq.ph$Q[same_markers], -->
<!--                         submaps_merged[[1]]$maps[[1]]$seq.ph$P[same_markers], -->
<!--                         submaps_merged[[1]]$maps[[1]]$seq.ph$Q[same_markers]) -->
<!-- ``` -->

<!-- The procedure produced a good map with the same phase configurations as those obtained using the sequential approach. Also, there was a significant reduction of the processing time: this approach produced the final map in approximately 46 seconds, while it took approximately 50 minutes to be finished using the sequential approach. -->

<!-- Even though the processing time advantage of this procedure, please be aware of the following points: -->

<!-- - As commented before, the approach is based on a known sequence order -->
<!-- - The time reduction is proportional to the number of batches and number of available cores -->
<!-- - The number of batches should be proportional to the density of the map, reliability and amount of information in each marker. Therefore, as the number of batches increase, the number of markers in each batch reduces, which may not be sufficient to identify the correct phase configuration. In general, a low number of batches presents more reliable results, as they will contain more information to build the submaps -->
<!-- - The same works for the definition of other parameters, such as rf and LOD thresholds. In general, you may inform lower values for these parameters when compared to the default HMM sequential approach  -->
<!-- - This procedure works well for a high number of markers (dense map); if you are dealing with a low number of markers, this approach may not produce good results. Even with good marker data, the amount of information in each marker for phasing reduces with the increase of distance. Then, the HMM will rely on wrong phase configurations during the propagation, which affects the re-estimation of recombination fractions and genotype probabilities. We strongly recommend that you perform as many tests as possible before continuing the analysis with your own data -->

## Map summary

After building two or more maps, one may want to compare summary statistics of those maps regarding the same chromosome, or even across chromosomes. A brief comparison can be done using the function `summary_maps`, which generates a table containing these statistics based on a list of `mappoly.map` objects:


```r
knitr::kable(summary_maps(MAPs))
```

# Genotype conditional probabilities

In order to use the genetic map in [QTLpoly](https://github.com/guilherme-pereira/QTLpoly), one needs to obtain the conditional probability of all possible 36 genotypes along the 12 linkage groups for all individuals in the full-sib population. Let us use the function `calc_genoprob_error`, which similarly to `est_full_hmm_with_global_error` allows the inclusion of a global genotyping error:


```r
genoprob.err <- vector("list", 12)
for(i in 1:12)
genoprob.err[[i]] <- calc_genoprob_error(input.map = MAP.err[[i]], error = 0.05)
```

Here, a global genotyping error of 5% was used. Each position of any object in the list `genoprob.err` contains two elements: an array of dimensions $36 \times number \; of \; markers \times  number \; of \; individuals$ and the position of the markers in the maps in centimorgans.  Let us display the results for all linkage groups in individual 1:


```r
ind <- 1
op <- par(mfrow = c(3, 4), pty = "s", mar=c(1,1,1,1)) 
for(i in 1:12)
{
  d <- genoprob.err[[i]]$map
  image(t(genoprob.err[[i]]$probs[,,ind]),
        col=RColorBrewer::brewer.pal(n=9 , name = "YlOrRd"),
        axes=FALSE,
        xlab = "Markers",
        ylab = "",
        main = paste("LG", i))
  axis(side = 1, at = d/max(d),
       labels =rep("", length(d)), las=2)
}
par(op)
```

In this figure, the x-axis represents the genetic map and the y-axis represents the 36 possible genotypes in the full-sib population. The color scale varies from dark purple (high probabilities) to light yellow (low probabilities). With the conditional probabilities computed, it is possible to use the object `genoprob.err` alongside with phenotypic data as the input of the software [QTLpoly](https://github.com/guilherme-pereira/QTLpoly), which is an under development software to map multiple QTLs in full-sib families of outcrossing autopolyploid species.

# Obtaining individual haplotypes

Once ready, the genotypic conditional probabilities can be used to recover any individual haplotype given the map (details described in @Mollinari2020). To generate this information, one may use the function `calc_homoprob` to account for the probabilities of each homologous, in all map positions for each individual. For example, let us view the homologous probabilities for chromosome 1 and individual 10:


```r
homoprobs = calc_homoprob(genoprob.err)
```


```r
plot(homoprobs, lg = 1, ind = 10)
```

Using this graphic, it is possible to identify regions of crossing-over occurrence, represented by the inversion of probability magnitudes between homologous from the same parent. It is also possible to view all chromosomes at the same time for any individual by setting the parameter `lg = "all"`. One may use this information to evaluate the quality of the map and repeat some processes with modifications, if necessary.

# Evaluating the meiotic process

MAPpoly also handles a function to evaluate the meiotic process that guided gamete formation on the studied population. Given genotype conditional probabilities, one may want to account for homologous pairing probabilities and detect the occurrence of preferential pairing, which is possible through the function `calc_prefpair_profiles`:


```r
prefpairs = calc_prefpair_profiles(genoprob.err)
```

The function returns an object of class `mappoly.prefpair.profiles`, which was saved as `prefpairs`. This object handles all information necessary to study pairing, such as the probability for each pairing configuration ($\psi$; see Mollinari and Garcia, 2019) inside each parent. For a more user-friendly visualization of the results, one may want to look at the `plot` output:


```r
plot(prefpairs)
#save.image("all.analysis.rda")
```

This graphic shows information about all pairing configurations and their probabilities, the proportion of bivalent/multivalent pairing, and also the p-value for preferential pairing test for all markers inside each parent.


# References

