---
title: "Copy number variation (CNV) analysis of hGIC samples"
output: html_notebook
author: "Gabriel Rosser"
---

We have previously observed a bias in the direction of the differentially methylated regions (DMRs) computed in individual patients, comparing the GIC line (typically $n=2$ different passages) with one of:

- iNSC
- fibroblast
- iAPC
- iOPC
 
The patients separate into two groups, one showing bias towards hypermethylation in GIC and vice versa. 

Here, we test whether this observation is related to CNVs in the genome. Hypothesis: deletions and amplifications are both associated with a reduced estimate of methylation. Therefore, those lines exhibiting hypomethylation have more major deletions or amplifications.

We approach this by estimating CNVs directly from methylation array data. We use the following packages:

- `minfi`: does all the heavy lifting, parsing `idat` files.
- `ChAMP`: nice wrapper over `minfi` that makes it easy to filter out bad probes. Might eliminate this if we can.
- `conumee`: estimates CNV in a single sample relative to a control (assumed to reflect an unperturbed genome).

For the control, we can use any one of the comparators listed above, or all of them. We'll try a few different settings.





Setup the data import directories. We simplify all this by using a metadata file generated in Python that lists all basenames, patient IDs and cell types.


```r
meta <- read.csv('methylation/hgic_metadata.csv', row.names = 1)

# for testing purposes: use fewer samples
meta <- meta[meta$patient_id %in% c(18, 19, 30),]

basenames <- file.path(data.dir.raid, 'methylation', rownames(meta))
```

Load the data. We're using `ChAMP` here to run some additional filtering, including removing probes with known SNPs and removing those below a minimum detection threshold.


```r
rgSet <- read.metharray(basenames, extended = T, force = T)
mSet <- preprocessSWAN(rgSet)

# We usually run ChAMP as follows:
# detP <- detectionP(rgSet)
# champLoad <- champ.filter(beta=NULL, M=mSet, detP = detP, pd=NULL, arraytype = 'EPIC')
# mSet <- champLoad$M
# However, in this context it introduces complications with the downstream conumee package. Therefore leave all probes in.
```

We now have a `MethylSet` object that has been preprocessed using the `SWAN` algorithm. This contains the estimated intensity signal for methylated and unmethylated versions of each probe, from which we can access the $\beta$ or $M$ values fairly straightforwardly.

Now we will use the package `conumee` to create an annotation, which divides the genome into bins and defines detailed regions of interest. We can use their pre-saved lists for this. We do need to filter this afterwards to remove any probes that are missing in our samples.



```r
data("detail_regions")  # supplied by conumee
data("exclude_regions")  # supplied by conumee
anno <- CNV.create_anno(array_type = 'EPIC', exclude_regions=exclude_regions, detail_regions = detail_regions)
```

```
## using genome annotations from UCSC
```

```
## getting EPIC annotations
```

```
##  - 844316 probes used
```

```
## importing regions to exclude from analysis
```

```
## importing regions for detailed analysis
```

```
## creating bins
```

```
##  - 53891 bins created
```

```
## merging bins
```

```
##  - 25735 bins remaining
```

```r
common_probes <- rownames(mSet)[is.element(rownames(mSet), names(anno@probes))]
anno@probes <- anno@probes[common_probes]
```

Finally, 'fit' our `MethylSet` to the annotation using `conumee` and run a series of queries, each of which compares a single GIC sample with a number of control samples.


```r
computeCNV <- function(query, control, anno, name=NULL) {
  x <- CNV.fit(query, control, anno, name = name)
  x <- CNV.bin(x)
  x <- CNV.detail(x)
  x <- CNV.segment(x)
}

minfi.data <- CNV.load(mSet)
```

```
## Warning in CNV.check(object): intensities are abnormally low (< 5000).
```

```r
cnv.res <- list()
for (pid in unique(meta$patient_id)) {
  gic.idx <- which((meta$patient_id == pid) & (meta$type == 'GBM'))
  ctrl.idx <- which((meta$patient_id == pid) & (meta$type != 'GBM'))
  cnv.res[[pid]] <- lapply(
    seq(1, length(gic.idx)), 
    function(x) computeCNV(minfi.data[gic.idx[x]], minfi.data[ctrl.idx], anno, name=sprintf("GIC%03d_%d", pid, x))
  )

}
```

# References

Feber A, Guilhamon P, Lechner M, et al. Using high-density DNA methylation arrays to profile copy number alterations. _Genome Biology_ 2014;15(2):R30. doi:10.1186/gb-2014-15-2-r30

[The `conumee` vignette](http://www.bioconductor.org/packages/devel/bioc/vignettes/conumee/inst/doc/conumee.html)

Maksimovic, J., Gordon, L., and Oshlack, A. (2012). SWAN: Subset-quantile Within Array Normalization for Illumina Infinium HumanMethylation450 BeadChips. _Genome Biology_ 13, R44.

