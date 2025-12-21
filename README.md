[![DOI](https://zenodo.org/badge/306393966.svg)](https://zenodo.org/badge/latestdoi/306393966)

# QDNAseq.hs1: QDNAseq bin annotation for hs.1

>QDNAseq bin annotation for the human genome build hs.1

This package provides QDNAseq bin annotations of size `1, 5, 10, 15, 30, 50, 100, 500 and 1000` kbp for the human genome build hs.1.The bin annotations are created using the steps mentioned in QDNAseq vignette and also [here](https://github.com/ccagc/QDNAseq/issues/59).


## Installation

Install the package from GitHub:

``` r
#Install the QDNAseq.hs1 package using remotes
remotes::install_github("GastroEsoLab/QDNAseq.hs1@main")
#or devtools
devtools::install_github("GastroEsoLab/QDNAseq.hs1@main")
```

## Use QDNAseq.hs1

``` r
library(QDNAseq)
library(QDNAseq.hs1)
bins <- getBinAnnotations(binSize=500, genome="hs1")
```

`QDNAseq.hs1` is adapted from [QDNAseq.hg38](10.5281/zenodo.4274555). Find more details about QDNAseq here: https://doi.org/doi:10.18129/B9.bioc.QDNAseq

`QDNAseq.hs1` has relatively strict mappability and filtering for artfactual regions. If you want more permissive annotations, please checkout [QDNAseq.chm13v2](https://github.com/RodrigoGM/QDNAseq.chm13v2).

## Make your own annotations
### Data & Tools Required
- 50mer hoffman mappability track from [UCSC](https://hgdownload.gi.ucsc.edu/gbdb/hs1/hoffmanMappability/k50.Umap.MultiTrackMappability.bw)
- Diploid samples from 1000 genomes for residual calculations (see bottom)
- Repeat Masker File from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz)

### 1. Create initial qDNASeq object
``` R
library(QDNAseq)
library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
library(rtracklayer)
library(GenomicRanges)
library(Biobase)
library(digest)
library(future.apply)

make_hs1_bins <- function(bin_size_kb, mappability_bw, bigwig_exec, genome, outdir) {
  
  message("\n====================================")
  message("Generating bins for: ", bin_size_kb, " kb")
  message("====================================")
  
  # Create bins
  bins <- createBins(bsgenome = genome, binSize = bin_size_kb)
  
  message("Calculating mappability...")
  
  bins$mappability <- calculateMappability(bins, bigWigFile = mappability_bw, bigWigAverageOverBed = bigwig_exec)
  
  # Define usable bins
  bins$use <- bins$chromosome %in% as.character(1:22) &
    bins$bases > 0 &
    bins$gc > 0 &
    bins$mappability > 0.5
  
  # Create annotated dataframe
  adf <- AnnotatedDataFrame(data = bins, varMetadata = data.frame(labelDescription = names(bins)))
  
  attr(adf, "QDNAseq") <- list(
    author = Sys.getenv("USER"),
    date = Sys.time(),
    organism = "Hsapiens",
    genome = "CHM13v2.0 / hs1",
    binsize = paste0(bin_size_kb, "kb"),
    mappability = "50mer 0 mismatch (GenMap)",
    url = NA,
    md5 = digest::digest(adf@data),
    sessionInfo = capture.output(sessionInfo())
  )
  
  outfile <- file.path(outdir, paste0("hs1.", bin_size_kb, "kb.SR50.rda"))
  save(adf, file = outfile, compress = "xz")
  
  message("Saved: ", outfile)
  return(outfile)
}

genome <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0
bin_sizes <- c(1000, 500, 100, 50, 30, 15, 10, 5, 1)
mappability_bw <- "/mnt/shera/bkw2118/make_QDNASeq.hs1/map/k50.Umap.MultiTrackMappability.bw"
bigwig_exec <- "/mnt/shera/bkw2118/make_QDNASeq.hs1/map/bigWigAverageOverBed"
outdir <- "/home/bkw2118/qDNASeq.hs1/"
dir.create(outdir, showWarnings = FALSE)

for (bin_size in bin_sizes) {
  make_hs1_bins(bin_size)
}

```

### 2. Calculate the Residuals against known diploid samples

``` R
# compute_CHM13_residuals.R
suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
  library(GenomicRanges)
})

options(future.globals.maxSize= 10 * 1024 ^ 3)

# -------- USER CONFIG ----------
bin_objects <- list.files('/home/bkw2118/qDNASeq.hs1/', pattern = '*.rda', full.names = T)   # change to the bin size you prefer
controls_dir <- "/mnt/shera/bkw2118/make_QDNASeq.hs1/bams/"              # folder with .bam remapped to hs1
bam_pattern <- "\\.bam$"                      # pattern to pick bam files
# --------------------------------

for (bins_rda in bin_objects){
  message("Processing bin object: ", bins_rda)
  out_residuals_rda <- gsub("\\.rda$", "_residuals.rda", bins_rda)
  load(bins_rda)

  # Get BAM file list
  bamfiles <- list.files(controls_dir, pattern=bam_pattern, full.names=TRUE)

  # Bin reads
  message("Binning Reads")
  readCounts <- binReadCounts(bins, bamfiles=bamfiles)
  
  # Filter bins
  readCountsFiltered <- applyFilters(readCounts, residual=FALSE, blacklist=FALSE)
  covs <- colSums(readCountsFiltered@assayData$counts)
  
  
  # We have enough samples so we can filter out low-coverage ones
  min_cov <- median(covs) * 0.25
  keep <- covs > min_cov
  if (sum(keep) < length(keep)) {
    readCountsFiltered <- readCountsFiltered[, keep]
  }

  # Calculate residuals compared to the expected CN State (all diploid)  
  message("Computing residuals with iterateResiduals()...")
  residvec <- iterateResiduals(readCountsFiltered)
  if (exists("adf")) {
    adf@data$residual <- residvec
    adf@data$use <- adf@data$use & !is.na(adf@data$residual)
    save(adf, file = out_residuals_rda, compress = "xz")
  } else {
    bins$residual <- residvec
    bins$use <- bins$use & !is.na(bins$residual)
    save(bins, file = out_residuals_rda, compress = "xz")
  }
  
  message("Residuals saved to: ", out_residuals_rda)
}
```

### 3. Create exclusion list based on repeats
``` R
# Load the repeat masker and merge with mappability data
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(MASS)
options(scipen = 999)

# Import repeat masker data
repeat_masker <- fread('/mnt/shera/bkw2118/make_QDNASeq.hs1/repeatMasker/hs1.repeatMasker.out')
colnames(repeat_masker) <- c("swScore", "milliDiv", "milliDel", "milliIns", "genoName", 
                             "genoStart", "genoEnd", "genoLeft", "strand", "repName", 
                             "repClass", "repStart", "repEnd", "repLeft", 
                             "id")
prob_repeats <- repeat_masker[grep('Satellite|Simple_repeat|Low_complexity|LINE|LTR|rRNA|tRNA|srpRNA', repeat_masker$repClass),]
repeat_GRanges <- GRanges(IRanges(prob_repeats$genoStart + 1, prob_repeats$genoEnd), 
                           seqnames = prob_repeats$genoName)
excl_regions <- reduce(repeat_GRanges)

# Filter out small regions
excl_regions <- excl_regions[width(excl_regions) > 1000]

# Convert into dataframe
excl_df <- data.frame(seqnames = as.character(seqnames(excl_regions)),
                      start = start(excl_regions),
                      end = end(excl_regions))
excl_df <- data.table(excl_df)

# Merge into the QDNASeq objects
qdna_files <- list.files('/home/bkw2118/QDNASeq.hs1', pattern = '*residuals.rda', full.names = TRUE)
for(file in qdna_files){
  load(file)
  bins <- adf@data
  
  binSize <- max(unique(bins$end - bins$start + 1))
  excl_df$bin_start <- floor((excl_df$start - 1) / binSize) * binSize + 1
  excl_df$bin_end <- excl_df$bin_start + binSize - 1
  excl_df$binName <- paste0(excl_df$seqnames, ":", excl_df$bin_start, "-", excl_df$bin_end)
  
  # Collapse by binName and get total number of bases excluded in each bin
  binned_excl <- excl_df[, .(sum(end - start, na.rm = T)/binSize), by = binName]
  binned_excl$exclude <- binned_excl$V1 >= 0.15
  qdna_binNames <- paste0('chr', bins$chromosome, ":", bins$start, "-", bins$end)
  
  bins$repetitive <- binned_excl$exclude[match(qdna_binNames, binned_excl$binName, nomatch = NA)]
  bins$repetitive[is.na(bins$repetitive)] <- FALSE
  bins$use <- (bins$mappability > 0.5) & (bins$gc > 30) & (bins$gc < 55) & (!bins$repetitive)
  adf@data <- bins
  
  ############# RENAME THE ADF TO THE FINAL NAME YOU WANT ##############
  hs1.5kb.SR50 <- adf
  ######################################################################
  save(hs1.5kb.SR50, file = gsub('residuals', 'final', file))
}
```

The Samples used for residual calculation were mostly derived from the 1000 genomes list

``` bash
HG01101
HG01204
HG01495
HG01522
NA06994
NA07051
NA11830
NA11832
NA11918
NA11919
NA11992
NA11994
NA12005
NA18489
NA18498
NA18504
NA18623
NA18624
NA18632
NA18633
NA18636
NA18867
NA18912
NA18960
NA18982
NA18984
NA18986
NA19058
NA19063
NA19064
NA19066
NA19116
NA19138
NA19474
NA19703
NA19707
NA19789
NA19901
```

