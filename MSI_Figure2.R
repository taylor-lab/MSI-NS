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

colo_tsg_dat <- fread("Colo_TSG.txt")
endo_tsg_dat <- fread("Endo_TSG.txt")

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


