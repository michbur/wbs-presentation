library(dplyr)
library(ggplot2)
library(reshape2)

load("./data/times2.RData")

times_dat <- do.call(rbind, lapply(times, function(i) 
  t(sapply(i, function(j) {
    c(slow = unname(j[["slow"]]["elapsed"]),
      quick = unname(j[["quick"]]["elapsed"]))
  }))
)) %>% 
  data.frame %>% 
  mutate(size = sort(rep(1L:5*10, 10))) %>% 
  melt(id.vars = "size") %>% 
  group_by(size, variable) %>% 
  summarise(value = mean(value)) %>% 
  ungroup %>% 
  mutate(variable = factor(variable, labels = c("Test permutacyjny", "QuiPT")))

cairo_ps("/home/michal/Dokumenty/overleaf/start/figures/QuiPT.eps", width = 2, height = 2)
ggplot(times_dat, aes(x = size, y = value, color = variable, shape = variable, linetype = variable)) +
  geom_point() +
  geom_line() +
  scale_y_continuous("Czas obliczeń [s]") +
  scale_x_continuous("Liczba testowanych n-gramów") +
  scale_color_manual("", values = c("red", "blue")) +
  scale_shape_discrete("") +
  scale_linetype_discrete("") +
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.34, 0.9),
        legend.direction = "vertical",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))
dev.off()

library(grid)
library(gridExtra)


best_enc <- AmyloGram::AmyloGram_model[["enc"]]
best_enc_aa <- do.call(rbind, lapply(1L:length(best_enc), function(i)
  data.frame(id = i, aa = toupper(best_enc[[i]]), stringsAsFactors = FALSE)
))
  
g_legend<-function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

gr_aa <- group_by(best_enc_aa, id) %>% 
  summarise(gr = paste0("{", paste0(aa, collapse = ""), "}")) 


ngram_freq <- read.csv("ngram_freq.csv")

ngram_freq_plot <- mutate(ngram_freq, decoded_name = gsub("_", "-", decoded_name)) %>%
  mutate(decoded_name = factor(decoded_name, levels = as.character(decoded_name)),
         amyloid = diff_freq > 0) %>%
  melt() %>%
  filter(variable %in% c("pos", "neg")) %>%
  droplevels %>%
  mutate(variable = factor(variable, labels = c("Amyloid", "Non-amyloid")))

for(i in 1L:6)
  levels(ngram_freq_plot[["decoded_name"]]) <- gsub(as.character(i), gr_aa[i, "gr"],
                                                    levels(ngram_freq_plot[["decoded_name"]]))

labels_colors <- c("black", "chartreuse3", "dodgerblue2", "firebrick1", "darkorange", "darkseagreen4", "darkorchid3")

gen_labels <- function(single_gr, x, gr_aa) {
  new_lab <- x
  for (other_gr in gr_aa[["gr"]][-single_gr])
    levels(new_lab) <- gsub(as.character(other_gr), paste0(rep(" ", nchar(other_gr)), collapse = ""), 
                            levels(new_lab), fixed = TRUE)
  levels(new_lab) <- gsub("-", " ", levels(new_lab), fixed = TRUE)
  new_lab
}

# generate list of labels where every element (amino acid group or dash) is present only once

all_labels <- c(eval({
  new_lab <- ngram_freq_plot[["decoded_name"]]
  levels(new_lab) <- gsub("[,A-Z\\{}]", " ", levels(new_lab))
  list(new_lab)
}),
lapply(1L:6, function(i) gen_labels(i, ngram_freq_plot[["decoded_name"]], gr_aa)))

levels(all_labels[[1]]) <- gsub("-", "X", levels(all_labels[[1]]))

# create series of plots where only one element (group of amino acids or dash) is plotted 

ngram_freq_plot2 <- filter(ngram_freq_plot, association != "Not found") %>% 
  droplevels %>% 
  mutate(variable = factor(variable, labels = c("Amyloid", "Nieamyloid")))

ngram_plots <- lapply(1L:7, function(i) {
  p <- ggplot(ngram_freq_plot2, aes(x = decoded_name, y = value)) +
    geom_bar(aes(fill = variable), position = "dodge", 
             stat = "identity", color = "black", size = 0.1,
             width = 0.7) +
    #geom_point(data = group_by(ngram_freq_plot2, decoded_name)  %>% filter(value == max(value)),
    #           aes(y = value + 0.004, shape = association), size = 2, stroke = 0.2, fill = "white") +
    scale_fill_manual("", values = c("white", "black")) +
    # scale_shape_manual("Experimentally validated motif:", breaks = c("Amyloidogenic", "Non-amyloidogenic"), 
    #                    values = c(21, 16, NA)) +
    scale_y_continuous("Częstość") +
    scale_x_discrete("", labels = all_labels[[i]]) + 
    my_theme +
    theme(axis.text.y = element_text(size = 6, colour = labels_colors[i], family = "mono", face = "bold"),
          legend.box = "vertical",
          panel.background = element_blank()) +
    coord_flip() 
})

amyl_legend <- g_legend(ngram_plots[[1]])

ngrams_plots_final <- lapply(1L:length(ngram_plots), function(i)
  if(i < 7) {
    arrangeGrob(ngram_plots[[i]] + 
                  theme(legend.position="none",
                        axis.title.x = element_text(size=8, vjust = -1, color = "white")),
                rectGrob(x = unit(0.5, "npc"), y = unit(0.5, "npc"), gp = gpar(col = "white")), 
                nrow = 2, heights = c(0.96, 0.04))
  } else {
    arrangeGrob(ngram_plots[[i]] + theme(legend.position="none"),
                amyl_legend, 
                nrow = 2, heights = c(0.96, 0.04))
    
  }
)

# combine plots
cairo_ps("/home/michal/Dokumenty/overleaf/start/figures/ngrams.eps", height = 3, width = 3)
for(i in 1L:7) {
  grid.draw(ngrams_plots_final[[i]])
}
dev.off()

source("script.R")  


cairo_ps("/home/michal/Dokumenty/overleaf/start/figures/SP2.eps", width = 3.2, height = 2.5)
rbind(mutate(pca_full[, 1L:3], alph = "Pełny alfabet"),
      mutate(pca_naive, alph = "Uproszczony alfabet")) %>% 
  ggplot(aes(x = PC1, y = PC2, color = taxon, shape = taxon, fill = taxon)) + 
  stat_density2d(aes(fill=taxon, alpha=..level..), 
                 color = "black", contour = TRUE, geom="polygon") +
  scale_linetype_discrete("") +
  scale_fill_manual("", values = c("orange", "blue"), labels = c("Inne eukarionty", "Plasmodiidae")) +
  scale_shape_discrete("") +
  scale_color_discrete("") +
  scale_alpha_continuous(range = c(0.25, 0.4)) +
  facet_wrap(~ alph, scales = "free", ncol = 2) +
  guides(alpha = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.725, 0.11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        legend.key.size = unit(0.12, "in"))
dev.off() 


ftir_res <- read.csv("ftir_amylogram_pasta2_literature.csv")

ftir_res_agg <- ftir_res %>% 
  group_by(amyl_db, amyl_ftir) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  mutate(amyl_db = factor(amyl_db, labels = c("Non-amyloid", "Amyloid")),
         amyl_ftir = factor(amyl_ftir, labels = c("Nieamyloid", "Amyloid"))) %>% 
  filter(amyl_db == "Non-amyloid",
         !is.na(amyl_ftir)) %>% 
  mutate(count = ifelse(count == 7, 4, count)) %>% 
  rbind(data.frame(amyl_db = "Non-amyloid", amyl_ftir = "Amyloid", count = 3)) %>% 
  mutate(frac = count/sum(count),
         count_nice = ifelse(count > 2, paste0(count, " peptydy"), paste0(count, " peptyd")),
         amyl_db = factor(c("Nieamyloid", 
                            "Amyloid", 
                            "Amyloid (literatura)")),
         amyl_db = factor(amyl_db, levels = levels(amyl_db)[c(3, 1, 2)]))

cairo_ps(filename = "/home/michal/Dokumenty/overleaf/start/figures/ftir.eps", width = 2.3, height = 2)
ggplot(ftir_res_agg, aes(x = amyl_ftir, y = frac, fill = amyl_db,
                         label = count_nice)) +
  geom_bar(stat = "identity") +
  geom_text(position = position_stack(vjust = 0.5), color = "black", size = 3) +
  scale_y_continuous("Frakcja peptydów") +
  scale_x_discrete("Wynik FTIR") +
  scale_fill_manual("", values = c(rgb(0, 156, 219, maxColorValue = 255),
                                   rgb(255, 65, 50, maxColorValue = 255),
                                   rgb(255, 170, 0, maxColorValue = 255))) +
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_blank(),
        plot.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        legend.key.size = unit(0.12, "in"))
dev.off()
