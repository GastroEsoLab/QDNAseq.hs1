# QDNAseq.hs1: QDNAseq bin annotation for hs.1

This package provides QDNAseq bin annotations of size `1, 5, 10, 15, 30, 50, 100, 500 and 1000` kbp for the human genome build hs.1.The bin annotations are created using the steps mentioned in QDNAseq vignette and also from the [hg38 QDNASeq](https://github.com/asntech/QDNAseq.hg38).

QDNAseq.hs1 is adapted from [QDNAseq.hg38](10.5281/zenodo.4274555). It was developed for use with [Songbird](https://github.com/GastroEsoLab/Songbird), so has relatively strict mappability and filtering for artfactual regions. If you want more permissive annotations, please checkout [QDNAseq.chm13v2](https://github.com/RodrigoGM/QDNAseq.chm13v2).

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

## Make your own annotations
In brief creating your own qDNASeq annotation object involves:
1. Extract mappability scores for *n*-mer sized reads and create the initial annotation object with gc content.
   - In our case we stick with the relatively conservative *50*mer window for mappability as Songbird uses fragments as short as 50 nucleotides for its true ploidy estimate.
2. Use known diploid samples to find regions in the genome with artificially higher or lower read density. We have gathered ~250 samples from the 1000 genomes project in order to calculate these residuals for read depth correction. This may or may not be useful depending on your sequencing platform and library preparation.
3. Create a *poor man's* encode exclusion list by finding difficult to align genomic regions and excluding them
   - In our case we flag satellite regions, retrotransposase regions, centromeric regions, tRNA, rRNA, srpRNA regions as difficult
   - Difficult regions which are more than 2000 nucleotides long are retained
4. Set the final `use` column. In this case, bins which do not include any difficult regions, and on average have a mappability score greater than 50, GC content greater than 30% and less than 70% are identified as useable. 

### Data & Tools Required
- hs1 reference genome from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz)
- 50mer hoffman mappability track from [UCSC](https://hgdownload.gi.ucsc.edu/gbdb/hs1/hoffmanMappability/k50.Umap.MultiTrackMappability.bw)
- Diploid samples from 1000 genomes for residual calculations (see bottom)
- Repeat Masker File from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz)
- the [bigWigAverageOverBed](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed) tool from UCSC (link points to linux binaries)

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

| | | | |
|---|---|---|---|
|ERR008834|SRR014155|SRR018029|ERR008835|
|SRR014156|SRR018030|ERR008838|SRR014157|
|SRR018031|ERR008839|SRR014159|SRR018032|
|ERR008841|SRR014160|SRR018033|ERR008843|
|SRR014161|SRR018034|ERR008941|SRR014162|
|SRR018035|ERR008942|SRR014163|SRR018036|
|ERR008957|SRR014164|SRR018103|ERR008958|
|SRR014165|SRR018104|ERR008959|SRR014175|
|SRR018105|ERR008960|SRR014176|SRR018106|
|ERR008963|SRR014177|SRR018107|ERR008964|
|SRR014178|SRR018108|ERR008972|SRR014179|
|SRR018109|ERR008974|SRR014180|SRR019043|
|ERR008975|SRR014181|SRR029726|ERR009041|
|SRR014182|SRR029727|ERR009051|SRR014183|
|SRR029728|SRR003380|SRR014184|SRR029729|
|SRR003381|SRR014185|SRR029730|SRR003382|
|SRR014186|SRR029731|SRR003383|SRR014187|
|SRR029732|SRR003384|SRR014188|SRR029733|
|SRR003385|SRR014189|SRR029825|SRR003386|
|SRR014190|SRR029826|SRR003428|SRR014191|
|SRR029827|SRR003429|SRR014192|SRR029828|
|SRR003430|SRR014193|SRR029829|SRR003431|
|SRR014194|SRR029830|SRR003432|SRR014195|
|SRR029831|SRR003433|SRR014196|SRR029832|
|SRR003442|SRR014197|SRR031307|SRR003443|
|SRR014198|SRR031308|SRR003444|SRR014199|
|SRR031338|SRR003445|SRR014200|SRR031339|
|SRR003446|SRR014201|SRR031340|SRR003447|
|SRR014202|SRR031341|SRR003448|SRR014203|
|SRR031342|SRR003660|SRR014204|SRR031343|
|SRR003661|SRR015466|SRR031344|SRR003662|
|SRR015467|SRR031345|SRR003663|SRR015468|
|SRR032377|SRR003664|SRR015469|SRR032378|
|SRR003665|SRR015470|SRR032379|SRR003666|
|SRR015471|SRR032380|SRR003927|SRR015472|
|SRR032381|SRR003928|SRR015473|SRR032382|
|SRR003929|SRR015983|SRR032383|SRR003930|
|SRR015984|SRR032384|SRR003931|SRR015985|
|SRR032387|SRR003932|SRR015987|SRR032388|
|SRR003933|SRR015988|SRR032389|SRR003934|
|SRR015989|SRR032390|SRR003935|SRR015990|
|SRR032391|SRR003936|SRR015991|SRR032392|
|SRR003937|SRR015992|SRR032393|SRR003938|
|SRR015993|SRR032394|SRR003939|SRR015994|
|SRR063405|SRR006215|SRR016228|SRR063406|
|SRR006216|SRR016229|SRR063407|SRR006217|
|SRR016230|SRR063408|SRR006218|SRR016231|
|SRR063409|SRR006273|SRR016232|SRR063410|
|SRR006274|SRR016233|SRR063411|SRR006275|
|SRR016234|SRR063412|SRR014127|SRR016235|
|SRR065212|SRR014128|SRR016405|SRR065216|
|SRR014129|SRR016406|SRR065218|SRR014130|
|SRR017498|SRR065219|SRR014131|SRR017499|
|SRR065222|SRR014132|SRR017500|SRR350130|
|SRR014133|SRR017501|SRR350133|SRR014134|
|SRR017502|SRR350141|SRR014135|SRR017503|
|SRR350144|SRR014136|SRR018023|SRR350145|
|SRR014137|SRR018024|SRR350149|SRR014138|
|SRR018025|SRR350152|SRR014152|SRR018026|
|SRR350162|SRR014153|SRR018027|SRR014154|
|SRR018028||||


