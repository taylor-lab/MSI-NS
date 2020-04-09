#!/usr/bin/env Rscript

pkgs = c('data.table',
         'HMMcopy',
         'ggplot2',
         'stringr',
         'parallel',
         'univOutl')

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

# load read counts (500kb windows)
# generated using readCounter from:
# https://github.com/shahcompbio/hmmcopy_utils
tfiles <- Sys.glob("hmmcopy/*.rg.md.23.wig")

# HMMcopy function
run_hmmcopy <- function(tfile) {
  
  mfile <- "map.500k.wig" # hmmcopy_utils mapCounter (b37)
  gfile <- "b37.gc.500k.wig" # hmmcopy_utils gcCounter (b37)
  
  # catverbose(tfile)
  tumour_copy <-
    correctReadcount(wigsToRangedData(tfile, gfile, mfile), verbose = F)
  
  default_param <- HMMsegment(tumour_copy, getparam = TRUE, verbose = F)
  
  new_param <- default_param
  new_param$strength <- 1e30
  
  # parameters tuned to focal heterzygous losses identified from WES
  new_param$mu <- c(-1.5, -0.75, 0, 0.375, 0.75, 1.5) 
  
  tumour_segments <- HMMsegment(tumour_copy, param = new_param, verbose = F)
  
  tumour_copy_dt <- suppressWarnings(as.data.table(tumour_copy))
  tumour_copy_dt[, state := tumour_segments$state]
  tumour_copy_dt[, space := factor(space, levels = c(1:22, "X"))]
  tumour_copy_dt <- tumour_copy_dt[!is.na(space)]
  tumour_copy_dt[, filename := tfile]
  copy(tumour_copy_dt)

}

##########################################################################################
##########################################################################################

sample_regex <- "(?<=[A-Z-]{2,8}[0-9]{2}T_)[A-Z0-9]{1,4}"
patient_regex <- "[A-Z-]{2,8}[0-9]{2}"

# run HMMcopy
hmmcopy_list <- lapply(tfiles, run_hmmcopy)
# hmmcopy_list <- parallel::mclapply(tfiles, run_hmmcopy, mc.cores = 2)

tumour_copy_dt <- rbindlist(hmmcopy_list)
tumour_copy_dt[, cid := stringr::str_extract(filename, "C[_-][A-Z0-9]{6}[_-]P[0-9]{3}[_-]d"), filename]
tumour_copy_dt[, cid := gsub("_", "-", cid), cid]

tumour_copy_dt <- merge(tumour_copy_dt, sample_name_key, by = "cid")

tumour_copy_dt[, sample := stringr::str_extract(pid, sample_regex), pid]
tumour_copy_dt[sample == "BULK", sample := "0"]
tumour_copy_dt[, sample := as.integer(sample)]
tumour_copy_dt[, patient := stringr::str_extract(pid, patient_regex), pid]

tumour_copy_dt[, state := state - 1]
tumour_copy_dt[, state := ifelse(state > 4, 4, state)] # cap amplifications

pts <- unique(tumour_copy_dt$patient)
# pts <-
#   c(
#     "CRC-MSI-01",
#     "CRC-MSI-02",
#     "CRC-MSI-03",
#     "CRC-MSI-04",
#     "CRC-MSI-05",
#     "CRC-MSI-06",
#     "EM01",
#     "EM03",
#     "EM04",
#     "EM05",
#     "EM06"
#   )

state_cols <- c(
  "0" = "#3182bd", # homozygous deletion
  "1" = "#9ecae1", # heterozygous loss
  "2" = NA, # diploid
  "3" = "#fcae91", # gain
  "4" = "#cb181d" # amplification
)

##########################################################################################
##########################################################################################

# DLRS QC
# out_qc <- c()

# univariate outlier detection
k_lsb = 2.5
lsb = LocScaleB(out_qc[sample != 0]$resid_sd, k = k_lsb, method = 'dq')
outliers <- out_qc[sample != 0][lsb$outliers]

# sample contamination
contamination <-
  data.table(
    patient = rep("EM06", 12),
    sample = c(0, 9, 17, 25, 33, 41, 49, 57, 65, 73, 81, 89)
  )

# flag cells without evidence of MSI
# normal flags from MSI_SC_Genotyping.R
tumour_copy_dt[, normal_flag := wes_sc_read_counts$normal_flag[match(cid, wes_sc_read_counts$cid)]]

# all_tcd_rle <- c()

# loop over patient IDs
for(pt in pts) {
  
  catverbose(pt)
  
  # remove probable normal cells, contaminated cells, and DLRS outliers
  tcd <- tumour_copy_dt[patient == pt & normal_flag == F] # normal flag
  tcd <- tcd[paste(patient, sample) %!in% contamination[, .(z=paste(patient, sample))]$z] # contamination
  tcd <- tcd[paste(patient, sample) %!in% outliers[, .(zz=paste(patient, sample))]$zz] # outliers
  
  # diagnostic plot
  # pdf(
  #   file = paste0("results/out/", pt, "_bias.pdf"),
  #   width = 10,
  #   height = 10
  # )
  # plotBias(tcd, pch = 20, cex = 0.5)
  # dev.off()
  
  dc <- dcast.data.table(tcd,
                         patient + space + start ~ pid,
                         value.var = "state")
  
  # hierarchical clustering for downstream visualization
  hc <- hclust(dist(t(dc[, -1:-3, with = F])), method = "ward.D2")
  tcd[, pid := factor(pid, levels = hc$labels[hc$order])]
  tcd <- tcd[order(pid)]
  sort_order <- tcd[, unique(sample)]
  sort_order <- c(0, setdiff(sort_order, 0))
  tcd[, sample := factor(sample, levels = sort_order)]
  
  # run length encoding
  tcd[, rle_id := rleid(sample, space, state)]
  tcd[, min_start := min(start), by = rle_id]
  tcd[, max_end := max(end), by = rle_id]
  tcd[, mean_copy := mean(copy, na.rm = T), by = rle_id]
  tcd[, state_sd := sd(copy, na.rm = T), by = list(sample, state)]
  tcd[, resid_sd := sd(copy - mean_copy, na.rm = T), by = list(sample)]
  tcd_rle = unique(tcd[, .(sample, space, state, rle_id, min_start, max_end, mean_copy, state_sd)])
  tcd_rle[, patient := pt]
  
  # out_qc <- rbind(out_qc, unique(tcd[, .(patient = pt, sample, resid_sd)]))
  # all_tcd_rle <- rbind(all_tcd_rle, tcd_rle)
  
  # }
  
  # plot contiguous segments
  p <- ggplot(tcd_rle[sample != "0" & space != "X"],
              aes(col = factor(state))) +
    geom_segment(data = tcd_rle[sample != "0" & space != "X"], aes(x = min_start, xend = max_end, y = 0, yend = 0), size = 2) +
    coord_cartesian(ylim = c(0, 0)) +
    facet_grid(sample ~ space, scales = "free_x",
               space = "free_x") +
    scale_color_manual("Copies",
                       values = state_cols,
                       labels = c("0", "1", "2", "3", "4", "4+")) +
    theme_classic() +
    theme(
      panel.spacing = unit(0, "lines"),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      # strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      legend.position = "top"
    ) +
    guides(color = guide_legend(nrow = 1))
  
  ggsave(
    p,
    filename = paste0(pt, "_hmmcopy_qc_pass_segs.pdf"),
    width = 20,
    height = 10
  )
  
}
