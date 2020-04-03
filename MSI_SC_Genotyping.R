#!/usr/bin/env Rscript

pkgs = c('data.table',
         'stringr',
         'plyr',
         'ggplot2',
         'RColorBrewer')

tmp <- lapply(pkgs, function (x) {
  suppressPackageStartupMessages(require(x, character.only = TRUE))
})
rm(tmp)

'%!in%' <- function(x, y)
  ! ('%in%'(x, y))
'%!like%' <- function(x, y)
  ! ('%like%'(x, y))

luq <- function(x) {return(length(unique(x)))}

catverbose <- function(...) {
  cat(format(Sys.time(), "%Y%m%d %H:%M:%S |"), ..., "\n")
}

##########################################################################################
##########################################################################################

# load Mutect2 force calls at mutant loci identified from WES
wes_sc_genotype_files <- Sys.glob("../*.out")
wes_sc_vcf <- lapply(wes_sc_genotype_files, function(x) fread(x, skip = "#CHROM"))
names(wes_sc_vcf) <- lapply(wes_sc_genotype_files, function(x) file_path_sans_ext(basename(x))[[1]])

wes_sc_genotypes = vector(mode = "list", length = length(wes_sc_vcf))
names(wes_sc_genotypes) <- names(wes_sc_vcf)

# parse VCF files
for (i in names(wes_sc_vcf)) {
  
  catverbose(i)
  
  tmp.dat <- wes_sc_vcf[[i]][ALT %!like% ","]
  tmp.out <- tmp.dat[, 1:9, with = F]
  cols <- grep("^C_", names(wes_sc_vcf[[i]]), value = T)
  
  for (j in cols) {
    tmp_col <- tstrsplit(tmp.dat[[j]], ":")[[2]]
    tmp_col_0 <- as.numeric(sapply(strsplit(tmp_col, ","), "[[", 1))
    tmp_col_1 <- as.numeric(sapply(strsplit(tmp_col, ","), "[[", 2))
    tmp.out[, paste0(j, "_", 0:1) := list(tmp_col_0, tmp_col_1)]
  }
  
  wes_sc_genotypes[[i]] <- tmp.out
  
}

wes_sc_read_counts = vector(mode = "list", length = length(wes_sc_vcf))
names(wes_sc_read_counts) <- names(wes_sc_vcf)

# loop over genotypes
# extract read counts
for (i in 1:length(wes_sc_genotypes)) {
  
  catverbose(names(wes_sc_genotypes)[i])
  
  sd.cols_0 <- grep("_0$", names(wes_sc_genotypes[[i]]), value = T)
  sd.cols_1 <- grep("_1$", names(wes_sc_genotypes[[i]]), value = T)
  
  tmp.out <-
    as.data.table(t(wes_sc_genotypes[[i]][, lapply(.SD, sum), .SDcols = sd.cols_0]), keep.rownames = T)
  tmp.out[, rn := gsub("_0$", "", rn)]
  setnames(tmp.out, "V1", "REF")
  tmp.out[, c("ALT") := t(wes_sc_genotypes[[i]][, lapply(.SD, sum), .SDcols = sd.cols_1])]
  tmp.out[, "Sample" := (names(wes_sc_genotypes)[i])]
  
  wes_sc_read_counts[[i]] <- tmp.out
  
}

wes_sc_read_counts <- rbindlist(wes_sc_read_counts)
wes_sc_read_counts <- wes_sc_read_counts[, .(Sample, cid = gsub("_", "-", rn), REF, ALT)]
wes_sc_read_counts[, TOTAL := REF + ALT]

# apply binomial test
# prior of 0.5 corresponds to random sampling of REF/ALT reads in diploid regions
wes_sc_read_counts[TOTAL > 0, pval := binom.test(
  x = ALT,
  n = TOTAL,
  p = 0.5,
  alternative = "less"
)$p.value, by = (wes_sc_read_counts[TOTAL > 0])]

# wes_sc_read_counts[, qval := p.adjust(pval, method = "bonferroni")]
# flag cells with > 10 total reads, < 5 MSI+ reads and p-value < 0.05
wes_sc_read_counts[, normal_flag := (TOTAL > 10 & ALT < 5 & pval < 0.05)]

# cells with < 10 total reads at mutant loci are indeterminate (likely due to insufficient tumor purity)
wes_sc_read_counts[TOTAL < 10, normal_flag := NA]

