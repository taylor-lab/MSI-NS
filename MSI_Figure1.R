#!/usr/bin/env Rscript

pkgs = c('data.table',
         'readxl',
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

# load sample data
msi.master <- as.data.table(readxl::read_excel("data/extended_data_tables.xlsx", sheet = 1))

# Fig 1a

# MSI score
tmp.dat = msi.master[, .(
  min = min(MSI_SCORE, na.rm = T),
  q1 = quantile(MSI_SCORE, probs = 0.25, na.rm = T),
  med = quantile(MSI_SCORE, probs = 0.5, na.rm = T),
  q3 = quantile(MSI_SCORE, probs = 0.75, na.rm = T),
  max = max(MSI_SCORE, na.rm = T)
), by = MSI]

p1 = ggplot(tmp.dat, aes(factor(MSI, levels = c("MSS", "UNK", "MSI")))) +
  geom_boxplot(aes(
    ymin = min,
    lower = q1,
    middle = med,
    upper = q3,
    ymax = max,
    fill = MSI
  ),
  stat = "identity") +
  scale_y_continuous(
    limits = c(0, 55),
    breaks = seq(0, 50, 10),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    "",
    values = c(
      "UNK" = "lightgrey",
      "MSS" = brewer.pal(n = 9, "Paired")[1],
      "MSI" = brewer.pal(n = 9, "Paired")[2]
    ),
    labels = rev(c("Indeterminate", "MSS", "MSI"))
  ) +
  ylab("MSI score") +
  xlab("") +
  theme_classic() +
  theme(
    aspect.ratio = 3,
    axis.line.x = element_blank(),
    axis.text = element_text(size = 15),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0, .25, 0, .25), 'in'),
    legend.position = "bottom"
  )

# Mutation rate
tmp.dat = msi.master[MSI != "UNK", .(
  min = log10(min(dmp_tmb + 0.01, na.rm = T)),
  q1 = log10(quantile(dmp_tmb + 0.01, probs = 0.25, na.rm = T)),
  med = log10(quantile(dmp_tmb + 0.01, probs = 0.5, na.rm = T)),
  q3 = log10(quantile(dmp_tmb + 0.01, probs = 0.75, na.rm = T)),
  max = log10(max(dmp_tmb + 0.01, na.rm = T))
), by = MSI]

p2 = ggplot(tmp.dat, aes(factor(MSI, levels = c("MSS", "MSI")))) +
  geom_boxplot(aes(
    ymin = min,
    lower = q1,
    middle = med,
    upper = q3,
    ymax = max,
    fill = MSI
  ),
  stat = "identity") +
  scale_y_continuous(
    limits = c(-2, 3),
    breaks = seq(-2, 3, 1),
    labels = c(0.01, 0.1, 1, 10, 100, 1000),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    "",
    values = c(
      "UNK" = "lightgrey",
      "MSS" = brewer.pal(n = 9, "Paired")[1],
      "MSI" = brewer.pal(n = 9, "Paired")[2]
    )
  ) +
  ylab("Mutation rate (mut/Mb)") +
  xlab("") +
  theme_classic() +
  theme(
    aspect.ratio = 4,
    axis.line.x = element_blank(),
    axis.text = element_text(size = 15),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0, .25, 0, .25), 'in'),
    legend.position = "none"
  )

# Indel ratio
tmp.dat = msi.master[MSI != "UNK" & (n_ins + n_del + n_snv > 9), .(
  min = min(indel_pct, na.rm = T),
  q1 = quantile(indel_pct, probs = 0.25, na.rm = T),
  med = quantile(indel_pct, probs = 0.5, na.rm = T),
  q3 = quantile(indel_pct, probs = 0.75, na.rm = T),
  max = max(indel_pct, na.rm = T)
), by = MSI]

p3 = ggplot(tmp.dat, aes(factor(MSI, levels = c("MSS", "MSI")))) +
  geom_boxplot(aes(
    ymin = min,
    lower = q1,
    middle = med,
    upper = q3,
    ymax = max,
    fill = MSI
  ),
  stat = "identity") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    "",
    values = c(
      "UNK" = "lightgrey",
      "MSS" = brewer.pal(n = 9, "Paired")[1],
      "MSI" = brewer.pal(n = 9, "Paired")[2]
    )
  ) +
  ylab("Indel ratio") +
  xlab("") +
  theme_classic() +
  theme(
    aspect.ratio = 4,
    axis.line.x = element_blank(),
    axis.text = element_text(size = 15),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0, .25, 0, .25), 'in'),
    legend.position = "none"
  )

# Fraction genome altered
tmp.dat = msi.master[MSI != "UNK", .(
  min = min(FGA, na.rm = T),
  q1 = quantile(FGA, probs = 0.25, na.rm = T),
  med = quantile(FGA, probs = 0.5, na.rm = T),
  q3 = quantile(FGA, probs = 0.75, na.rm = T),
  max = max(FGA, na.rm = T)
), by = MSI]

p4 = ggplot(tmp.dat, aes(factor(MSI, levels = c("MSS", "MSI")))) +
  geom_boxplot(aes(
    ymin = min,
    lower = q1,
    middle = med,
    upper = q3,
    ymax = max,
    fill = MSI
  ),
  stat = "identity") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    "",
    values = c(
      "UNK" = "lightgrey",
      "MSS" = brewer.pal(n = 9, "Paired")[1],
      "MSI" = brewer.pal(n = 9, "Paired")[2]
    )
  ) +
  ylab("Fraction genome altered") +
  xlab("") +
  theme_classic() +
  theme(
    aspect.ratio = 4,
    axis.line.x = element_blank(),
    axis.text = element_text(size = 15),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0, .25, 0, .25), 'in'),
    legend.position = "none"
  )

p1 + p2 + p3 + p4 + plot_layout(ncol = 4)

# Fig 1b

ggplot(msi.master[MSI == "MSI"],
       aes(x = FGA,
           fill = factor(FGA_bin))) +
  scale_fill_manual(
    "Copy number state",
    values = c(
      "Diploid" = "lightgrey",
      "Near Diploid" = sequential_hcl(n = 5, palette = "Reds 2")[3],
      "Moderate CNA" = sequential_hcl(n = 5, palette = "Reds 3")[2],
      "Extensive CNA" = sequential_hcl(n = 5, palette = "Reds")[1]
    )
  ) +
  geom_histogram(
    color = 'black',
    bins = 100,
    binwidth = 0.01,
    boundary = 0,
    closed = "left"
  ) +
  xlab("Fraction of genome altered") +
  ylab("Number of MSI cases") +
  scale_x_continuous(
    breaks = seq(0, 1, 0.1),
    labels = seq(0, 1, 0.1),
    expand = c(0.005, 0.005)
  ) +
  scale_y_continuous(limits = c(0, 150),
                     expand = c(0, 0)) +
  theme_classic() +
  theme(
    aspect.ratio = 0.5,
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0, .25, 0, .25), 'in'),
    legend.position = "none"
  )

# Fig 1c

fga_breakdown <-
  msi.master[MSI == "MSI", .N,
             by = list(FGA_bin, MSI_class)][order(MSI_class, FGA_bin)][, .(FGA_bin, N, pct = N / sum(N)),
                                                                       by = MSI_class %!like% 'Other']
fga_breakdown[, total := sum(N), by = MSI_class]
fga_breakdown[, MSI_class := ifelse(MSI_class, "MSI-assoc", "Other")]

ggplot(fga_breakdown,
       aes(
         x = (factor(
           MSI_class, levels = c("Other", "MSI-assoc")
         )),
         y = pct,
         fill = factor(FGA_bin, levels = rev(
           c("Diploid",
             "Near Diploid",
             "Moderate CNA",
             "Extensive CNA")
         ))
       )) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_y_continuous(
    expand = c(0.01, 0.01),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_fill_manual(
    "Copy number state",
    values = c(
      "Diploid" = "lightgrey",
      "Near Diploid" = sequential_hcl(n = 5, palette = "Reds 2")[3],
      "Moderate CNA" = sequential_hcl(n = 5, palette = "Reds 3")[2],
      "Extensive CNA" = sequential_hcl(n = 5, palette = "Reds")[1]
    ),
    labels = c("Diploid",
               "Near Diploid",
               "Moderate CNA",
               "Extensive CNA")
  ) +
  xlab("") +
  ylab("% MSI samples") +
  theme_classic() +
  theme(aspect.ratio = 2) +
  theme(
    plot.title = element_text(size = 25),
    axis.text = element_text(size = 15),
    axis.text.x = element_text(
      size = 15,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0, .25, 0, .25), 'in'),
    legend.position = "right"
  ) +
  guides(fill = guide_legend(override.aes = list(
    fill = c(
      "lightgrey",
      sequential_hcl(n = 5, palette = "Reds 2")[3],
      sequential_hcl(n = 5, palette = "Reds 3")[2],
      sequential_hcl(n = 5, palette = "Reds")[1]
    )
  )))

# Fig 1d

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
                              by = list(Cancer_Type, MSI_class %!like% 'Other', MSI)][order(Cancer_Type)]

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
    labels = c("MSI-assoc cancers",
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

# Fig 1e

# WGD rates by MSI class
ggplot(msi.master[MSI != "UNK", .N,
                  list(MSI_class %like% 'Other', MSI, WGD)][, .(WGD,
                                                               N,
                                                               label = paste(MSI_class, MSI),
                                                               pct = N / sum(N)), by = list(MSI_class, MSI)],
       aes(
         x = factor(
           label,
           levels = c("FALSE MSS", "FALSE MSI", "TRUE MSS", "TRUE MSI"),
           labels = c("MSI-assoc MSS", "MSI-assoc MSI", "Other MSS", "Other MSI")
         ),
         y = 100 * pct,
         fill = factor(WGD, levels = rev(c("TRUE", "FALSE")))
       )) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual("WGD",
                    values = c("FALSE" = "lightgrey",
                               "TRUE" = "black")) +
  scale_y_continuous(
    limits = c(0, 50),
    expand = c(0, 0),
    breaks = seq(0, 50, by = 10),
    labels = seq(0, 50, by = 10)
  ) +
  xlab("") +
  ylab("% of cases with WGD") +
  theme_classic() +
  theme(
    aspect.ratio = 2,
    axis.text = element_text(size = 20),
    axis.text.x = element_text(
      size = 20,
      vjust = 0.5,
      hjust = 1,
      angle = 90
    ),
    axis.title = element_text(size = 25),
    axis.ticks.length = unit(.25, "cm"),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.grid = element_blank(),
    legend.position = "none"
  )


