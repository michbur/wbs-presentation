library(dplyr)
library(ggplot2)
library(gridExtra)
source("my_ggplot_theme.R")


load("fr_benchmark_full.RData")
load("pca_results.RData")

plot_pca <- function(pca)
  ggplot(pca, aes(x = PC1, y = PC2, color = taxon, shape = taxon, fill = taxon)) + 
  stat_density2d(aes(fill=taxon, alpha=..level..), 
                 color = "black", contour = TRUE, geom="polygon") +
  scale_linetype_discrete("") +
  scale_fill_manual("", values = c("orange", "blue"), labels = c("Inne eukarionty", "Plasmodiidae")) +
  scale_shape_discrete("") +
  scale_color_discrete("") +
  scale_alpha_continuous(range = c(0.25, 0.4)) +
  guides(alpha = FALSE) +
  my_theme

dat <- unlist(fr_benchmark, recursive = FALSE) %>% 
  unlist(recursive = FALSE) %>% 
  do.call(rbind, .) %>% 
  lapply(unlist) %>% 
  data.frame() %>% 
  rename(MCC = MCC.MCC, AUC = AUC.AUC) %>% 
  mutate(classifier = "signalP",
         id = rep(sort(rep(1L:20, 2)), 9),
         sp = factor(sp, labels = c(20, 2, 4)),
         cs = factor(cs, labels = c(20, 2, 4))) %>% 
  mutate(sp = factor(sp, levels = levels(sp)[c(2, 3, 1)]),
         cs = factor(cs, levels = levels(cs)[c(2, 3, 1)]))

cmp_dat <- cbind(filter(dat, !plas),
                 filter(dat, plas) %>% 
                   select(MCC, AUC, mscsd) %>% 
                   rename(MCC_p = MCC, AUC_p = AUC, mscsd_p = mscsd)) %>% 
  filter(sp != 4, cs != 4)

p1 <- plot_pca(pca_full) +
  ggtitle("A")

p2 <- plot_pca(pca_naive) +
  ggtitle("B")

p3 <- ggplot(cmp_dat, aes(x = MCC, y = MCC_p)) +
  geom_point(aes(fill = sp), size = 2, shape = 21, color = "white") +
  geom_point(aes(color = cs), size = 3, shape = 1, stroke = 1) +
  my_theme +
  scale_x_continuous("MCC (other eukaryotes)") +
  scale_y_continuous("MCC (Plasmodiidae)") +
  scale_fill_manual("Length of the alphabet\n(signal peptide)", values = c("dodgerblue", "firebrick1")) +
  scale_color_manual("Length of the alphabet\n(cleavage site)", values = c("dodgerblue", "firebrick1")) +
  ggtitle("C")

cairo_ps("./figures/SP.eps", height = 5.4, width = 6.5)
grid.newpage()
grid.draw(arrangeGrob(arrangeGrob(p1, p2, ncol = 2), p3, ncol = 1, heights = c(0.45, 0.55)))
dev.off()
