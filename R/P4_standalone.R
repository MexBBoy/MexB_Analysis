# MexB pocket volume vs how deep the bound ligand sits.
# Self-contained: paste the whole file into R and run. Needs ggplot2;
# ggrepel is optional and only improves label placement.

library(ggplot2)

d <- data.frame(
  ligand = c("ampicillin", "DDM x3", "CYMAL-7", "chloramphenicol",
             "DDM (2V50)", "DDM (3W9I)", "EPI", "LMNG"),
  pdb    = c("this work", "this work", "21FO", "21FP",
             "2V50", "3W9I", "3W9J", "6IIA"),
  depth  = c(62.80, 52.38, 33.53, 62.13, 61.63, 62.13, 58.69, 52.80),
  volume = c(1850,  2068,  2037,  2132,  1953,  2209,  1868,  2008),
  atoms  = c(25,    105,   36,    20,    35,    35,    49,    69),
  ours   = c(TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

# poster colours, darkened at constant hue to stay legible on white
pal <- c("ampicillin"      = "#078A08",   # poster #1EFF21
         "DDM x3"          = "#CF13CF",   # poster #FF29FF  (DDM #1)
         "DDM (2V50)"      = "#986598",   # poster #FFB0FF  (DDM #2)
         "DDM (3W9I)"      = "#C24A8B",   # poster #FF6DBC  (legend DDM)
         "LMNG"            = "#2F54FF",   # poster #3E61FF
         "CYMAL-7"         = "#A1685E",   # poster #FFAB9C
         "EPI"             = "#767676",   # poster #A3A3A3
         "chloramphenicol" = "#000000")   # poster black

cat(sprintf("volume vs depth       r = %+.2f\n", cor(d$depth, d$volume)))
cat(sprintf("volume vs ligand size r = %+.2f\n", cor(d$atoms, d$volume)))

p <- ggplot(d, aes(depth, volume, colour = ligand)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "#104862", linetype = "22", linewidth = 0.6) +
  geom_point(aes(size = atoms, shape = ours)) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_size_area(max_size = 9, guide = "none") +   # area tracks ligand size
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  scale_x_continuous(limits = c(25, 80), breaks = seq(30, 70, 10)) +
  labs(x = expression("depth into the pocket from the periplasmic entrance ("
                      * ring(A) * ")"),
       y = expression("ligand-free pocket volume (" * ring(A)^3 * ")")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.line          = element_line(colour = "#9fb0b8"))

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(aes(label = ligand), size = 3.6, seed = 1,
                                    box.padding = 0.5, show.legend = FALSE)
} else {
  p <- p + geom_text(aes(label = ligand), size = 3.6, vjust = -1.4,
                     show.legend = FALSE)
}

print(p)
