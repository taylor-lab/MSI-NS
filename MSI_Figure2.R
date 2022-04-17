#!/usr/bin/env Rscript

pkgs = c('data.table',
         'stringr',
         'ggplot2',
         'RColorBrewer',
         'colorspace',
         'patchwork',
         'ggbeeswarm')

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

colo_tsg_dat <- fread("data/Colo_TSG.txt")
endo_tsg_dat <- fread("data/Endo_TSG.txt")

# Colorectal TSG biallelic inactivation

gg1 <- ggplot(colo_tsg_dat,
       aes(
         x = lab,
         y = pct,
         fill = factor(mode, levels = rev(c("HIT2", "CNLOH", "LOH", "DEL", "WT intact")))
       )) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_y_continuous(
    expand = c(0.0, 0.0),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  facet_wrap(~ factor(MSI, levels = c("MSS", "MSI")), scales = "free_x") +
  scale_fill_manual(
    "Mechanism",
    values = rev(
      c(
        "WT intact" = "white",
        "HIT2" = "dimgrey",
        "CNLOH" = sequential_hcl(n = 5, palette = "Oslo")[2],
        "LOH" = sequential_hcl(n = 5, palette = "Oslo")[3],
        "DEL" = sequential_hcl(n = 5, palette = "Oslo")[4]
      )
    ),
    breaks = rev(c("HIT2", "CNLOH", "LOH", "DEL")),
    labels = rev(c("2nd hit", "CNLOH", "LOH", "HomDel"))
  ) +
  xlab("") +
  ylab("% Biallelic inactivation") +
  theme_classic() +
  theme(aspect.ratio = 2) +
  theme(
    plot.title = element_text(size = 25),
    plot.subtitle = element_text(size = 15, face = "italic"),
    axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    plot.margin = unit(c(0, .25, 0, .25), 'in'),
    legend.position = "right"
  )

# Endometrial TSG biallelic inactivation

gg2 <- ggplot(endo_tsg_dat,
       aes(
         x = lab,
         y = pct,
         fill = factor(mode, levels = rev(c("HIT2", "CNLOH", "LOH", "DEL", "WT intact")))
       )) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_y_continuous(
    expand = c(0.0, 0.0),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  facet_wrap(~ factor(MSI, levels = c("MSS", "MSI")), scales = "free_x") +
  scale_fill_manual(
    "Mechanism",
    values = rev(
      c(
        "WT intact" = "white",
        "HIT2" = "dimgrey",
        "CNLOH" = sequential_hcl(n = 5, palette = "Oslo")[2],
        "LOH" = sequential_hcl(n = 5, palette = "Oslo")[3],
        "DEL" = sequential_hcl(n = 5, palette = "Oslo")[4]
      )
    ),
    breaks = rev(c("HIT2", "CNLOH", "LOH", "DEL")),
    labels = rev(c("2nd hit", "CNLOH", "LOH", "HomDel"))
  ) +
  xlab("") +
  ylab("% Biallelic inactivation") +
  theme_classic() +
  theme(aspect.ratio = 2) +
  theme(
    plot.title = element_text(size = 25),
    plot.subtitle = element_text(size = 15, face = "italic"),
    axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    plot.margin = unit(c(0, .25, 0, .25), 'in'),
    legend.position = "right"
  )

gg1 + gg2 + plot_layout(ncol = 2)

# Fig 2b

msi_cnloh_pcts <- fread("data/MSI_CNLOH_pcts.txt") # Percent of TSG biallelic inactivation due to CNLOH

ggplot(msi_cnloh_pcts[mode == "CNLOH"],
       aes(
         x = factor(
           MSI_class,
           levels = c("MSI-assoc MSS", "MSI-assoc MSI", "Non-MSI-assoc MSS", "Non-MSI-assoc MSI")
         ),
         y = 100 * pct
       )) +
  geom_bar(stat = "identity", width = 0.8, fill = "black") +
  geom_errorbar(aes(ymin = 100 * pct, ymax = 100 * pct.hi), width = 0.5, color = "black") +
  geom_errorbar(aes(ymin = 100 * pct.lo, ymax = 100 * pct), width = 0.5, color = "white") +
  scale_fill_manual("CNLOH",
                    values = c("CNLOH" = "black")) +
  scale_y_continuous(
    limits = c(0, 50),
    expand = c(0, 0),
    breaks = seq(0, 50, by = 10),
    labels = seq(0, 50, by = 10)
  ) +
  xlab("") +
  ylab("% of CNA-driven biaallelic\ninactivation due to CN-LOH") +
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

# Fig 2c

msi_loss_pcts <- fread("data/MSI_loss_pcts.txt") # Percent of CNAs due to genomic loss

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

sc_pct_fga <- fread("data/SC_pct_FGA.txt") # Percent FGA per single cell
pct_inc_clone <- fread("data/SC_pct_inc_clone.txt") # Percent single cells indistinguishable from incumbent clone 

p1 <- ggplot(sc_pct_fga,
             aes(x = patient, y = 100 * fga)) +
  geom_quasirandom(
    size = 3,
    alpha = 0.8,
    color = 'lightblue',
    width = 0.25
  ) +
  geom_violin(lwd = 1.1) +
  ylim(c(0, 80)) +
  xlab("") +
  ylab("Perc autosomal\n genome altered\n per single cell") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.ticks.length = unit(.25, "cm"),
    panel.grid = element_blank(),
  )

# FGA marginal histogram
fga_hist = hist(sc_pct_fga$fga,
                n = 40,
                xlim = c(0, 0.8))

barplot(
  fga_hist$counts,
  col = "black",
  border = "white",
  xlab = "Freq",
  axis.lty = 1,
  horiz = T,
)

p2 <- ggplot(pct_inc_clone,
             aes(x = patient, y = 100 * pct)) +
  geom_bar(stat = "identity",
           width = 0.9,
           fill = "black") +
  xlab("Patients (n=11)") +
  ylab("Perc\n incumbent\n clone") +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    panel.grid = element_blank(),
  ) +
  scale_y_reverse(limits = c(100, 0), 
                  breaks = c(0, 50, 100))

p1 / p2 + plot_layout(heights = c(2, 1))


# Fig 2e

sc_pt_crc05 <- fread("data/SC_CRC05.txt") # processed single-cell copy number profiles for patient CRC05

state_cols <- c(
  "0" = "#3182bd", # homozygous deletion
  "1" = "#9ecae1", # heterozygous loss
  "2" = "white", # diploid
  "3" = "#fcae91", # gain
  "4" = "#cb181d" # amplification
)

ggplot(sc_pt_crc05,
       aes(col = factor(state))) +
  geom_segment(data = sc_pt_crc05,
               aes(
                 x = min_start,
                 xend = max_end,
                 y = 0,
                 yend = 0
               ),
               size = 2) +
  coord_cartesian(ylim = c(0, 0)) +
  facet_grid(sample ~ chr, scales = "free_x",
             space = "free_x") +
  scale_color_manual("Copies",
                     values = state_cols,
                     labels = c("0", "1", "2", "3", "4+")) +
  theme_classic() +
  ylab("Individual tumor cells") +
  xlab("Chromosomes") +
  theme(
    panel.spacing = unit(0, "lines"),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "bottom",
  ) +
  guides(color = guide_legend(nrow = 1))




