#!/usr/bin/env Rscript

pkgs = c('data.table',
         'stringr',
         'ggplot2',
         'RColorBrewer',
         'colorspace',
         'patchwork')

tmp <- lapply(pkgs, function (x) {
  suppressPackageStartupMessages(require(x, character.only = TRUE))
})
rm(tmp)

'%!in%' <- function(x, y)
  ! ('%in%'(x, y))
'%!like%' <- function(x, y)
  ! ('%like%'(x, y))

##########################################################################################
##########################################################################################

# Fig 2a

msi_quantile_dat = msi.master[MSI != "UNK",
                              .(
                                hap = quantile(HAP, probs = 0.5, na.rm = T),
                                hap.lo = quantile(HAP, probs = 0.25, na.rm = T),
                                hap.hi = quantile(HAP, probs = 0.75, na.rm = T),
                                msi = quantile(MSI_SCORE, probs = 0.5, na.rm = T),
                                msi.lo = quantile(MSI_SCORE, probs = 0.25, na.rm = T),
                                msi.hi = quantile(MSI_SCORE, probs = 0.75, na.rm = T),
                                purity = quantile(Purity, probs = 0.5, na.rm = T),
                                purity.lo = quantile(Purity, probs = 0.25, na.rm = T),
                                purity.hi = quantile(Purity, probs = 0.75, na.rm = T),
                                N = luq(Tumor_Sample_Barcode)
                              ),
                              by = list(Cancer_Type, MSI_class %!like% 'Non-Lynch', MSI)][order(Cancer_Type)]

ggplot(msi_quantile_dat[N > 10 & MSI == "MSI"], aes(x = msi,
                                                    y = 100 * hap)) +
  geom_point(size = 5,
             alpha = 0.8,
             pch = 21,
             aes(fill = factor(MSI_class,
                               levels = c("TRUE", "FALSE")))) +
  geom_segment(
    data = test[N > 10 & MSI == "MSI"],
    color = "black",
    aes(
      y = 100 * hap,
      yend = 100 * hap,
      x = msi.lo,
      xend = msi.hi
    )
  ) +
  geom_segment(
    data = test[N > 10 & MSI == "MSI"],
    color = "black",
    aes(
      x = msi,
      xend = msi,
      y = 100 * hap.lo,
      yend = 100 * hap.hi
    )
  ) +
  scale_fill_manual(
    "",
    values = c("TRUE" = "#FDE330",
               "FALSE" = 'dimgrey'),
    labels = c("LS-associated tumors",
               "Other")
  ) +
  scale_radius(
    "Count",
    range = c(4, 16),
    trans = "sqrt",
    breaks = c(10, 50, 150)
  ) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  ylab("% haploid genome") +
  xlab("MSIsensor score") +
  scale_x_continuous(limits = c(12, 45), 
                     breaks = seq(15, 45, by = 10)) +
  ylim(c(0, 20)) +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 25),
    axis.ticks.length = unit(.25, "cm"),
    legend.text = element_text(size = 15),
    legend.position = c(0.9, 0.9)
  )

# Fig 2c

msi_loss_pcts <- fread("data/MSI_loss_pcts.txt")

ggplot(msi_loss_pcts[],
       aes(x = MSI_class, fill = clonality)) +
  geom_pointrange(
    position = position_dodge(width = 0.5),
    size = 1,
    shape = 21,
    color = "black",
    aes(y = pct,
        ymin = pct.lo,
        ymax = pct.hi)
  ) +
  scale_fill_manual("",
                    values = c("Subclonal" = brewer.pal(n = 9, "Paired")[1],
                               "Clonal" = brewer.pal(n = 9, "Paired")[2])
  ) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     breaks = seq(0, 10, by = 0.25),
                     labels = seq(0, 10, by = 0.25)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = 'red') +
  xlab("") +
  ylab("% CNAs due to genomic loss") +
  theme_classic() +
  theme(
    aspect.ratio = 2.5,
    axis.text.x = element_text(size = 20, hjust = 1, vjust = 0.5, angle = 90),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 25),
    axis.ticks.length = unit(.25, "cm"),
    legend.text = element_text(size = 15),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# Fig 2d

dndscv_dat <- fread("data/MSI_dndscv.txt")

ggplot(dndscv_dat,
       aes(x = zygosity_class, 
           fill = name)) +
  geom_pointrange(
    position = position_dodge(width = 0.5),
    size = 2,
    shape = 21,
    color = "black",
    aes(y = mle,
        ymin = cilow,
        ymax = cihigh)
  ) +
  facet_wrap(~ Ploidy, scales = "free_x") +
  scale_fill_manual("",
                    values = c("wmis" = "lightgrey",
                               "wtru" = "dimgrey"),
                    labels = c("Missense", "Truncating")
  ) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.5),
                     breaks = seq(0, 1.5, by = 0.25),
                     labels = seq(0, 1.5, by = 0.25)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("") +
  ylab("dN/dS") +
  theme_classic() +
  theme(
    aspect.ratio = 4,
    axis.text.x = element_text(size = 20, hjust = 1, vjust = 0.5, angle = 90),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 15),
    axis.ticks.length = unit(.25, "cm"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

