# ---------------------------------------------------------------------------
# P4: MexB pocket volume against how deep the bound ligand sits
#
# Reproduces results/figures/poster/P4_ligand_size_vs_pocket.png in ggplot2.
#
#   Rscript R/P4_pocket_depth.R
#
# Input : results/tables/published_pockets.csv   (written by
#         scripts/published_pockets.py; committed, so this runs standalone)
# Output: results/figures/poster/P4_R.pdf and .png
#
# Needs: ggplot2. ggrepel is used for the labels if installed; without it the
# script falls back to fixed label offsets and still runs.
# ---------------------------------------------------------------------------

library(ggplot2)

# find the repo root by walking up from the working directory, so the script
# runs from anywhere (repo root, R/, or an IDE session)
repo <- getwd()
for (i in 1:4) {
  if (file.exists(file.path(repo, "results", "tables",
                            "published_pockets.csv"))) break
  repo <- dirname(repo)
}
csv <- file.path(repo, "results", "tables", "published_pockets.csv")
if (!file.exists(csv)) stop("published_pockets.csv not found - run this from ",
                            "the MexB_Analysis repo")
d   <- read.csv(csv, stringsAsFactors = FALSE)

# ---- what goes on the plot -------------------------------------------------
# Ligand-bound structures only. The single apo entry (6T7S) is dropped: it is
# 4.5 A and much smaller, and scoring it as "0 heavy atoms" manufactures a
# correlation out of one low-resolution point.
d <- subset(d, ligand_heavy_atoms > 0 & !is.na(depth_from_entrance_A))

# Volume is quoted at a 16 A sphere. Columns for 14, 18 and 20 A are in the
# same file - a single radius is not meaningful on its own, so check that the
# conclusion survives the others before quoting any one of them.
d$volume <- d$volume_r16_A3
d$depth  <- d$depth_from_entrance_A

# ---- substrate names and the poster's colours ------------------------------
# Sampled from the conserved-residues legend and table header of
# Combio_Poster_20260828.pdf, then darkened at constant hue to clear 4.5:1 on
# white (the poster sets them on dark panels, where they run 1.4-2.5:1).
# 3W9I is not on the poster; it takes the pink the legend uses for DDM.
lig <- c(Amp_MexB_20260826   = "ampicillin",
         MexB_DDM_3_20260730 = "DDM x3",
         `2V50`              = "DDM (2V50)",
         `3W9I`              = "DDM (3W9I)",
         `6IIA`              = "LMNG",
         `21FO`              = "CYMAL-7",
         `3W9J`              = "EPI",
         `21FP`              = "chloramphenicol")

pal <- c(`ampicillin`      = "#078A08",   # poster #1EFF21
         `DDM x3`          = "#CF13CF",   # poster #FF29FF  (DDM #1)
         `DDM (2V50)`      = "#986598",   # poster #FFB0FF  (DDM #2)
         `DDM (3W9I)`      = "#C24A8B",   # poster #FF6DBC  (legend DDM)
         `LMNG`            = "#2F54FF",   # poster #3E61FF
         `CYMAL-7`         = "#A1685E",   # poster #FFAB9C
         `EPI`             = "#767676",   # poster #A3A3A3
         `chloramphenicol` = "#000000")   # poster black

d$ligand <- factor(lig[d$pdb], levels = names(pal))

# ---- the number in the callout ---------------------------------------------
r_depth <- cor(d$depth, d$volume)
r_size  <- cor(d$ligand_heavy_atoms, d$volume)
cat(sprintf("volume vs depth      r = %+.3f  (n = %d)\n", r_depth, nrow(d)))
cat(sprintf("volume vs ligand size r = %+.3f\n", r_size))

# ---- plot ------------------------------------------------------------------
p <- ggplot(d, aes(depth, volume, colour = ligand)) +
  geom_smooth(method = "lm", se = FALSE, colour = "#104862",
              linetype = "22", linewidth = 0.6, alpha = 0.55,
              formula = y ~ x) +
  geom_point(size = 4) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_continuous(limits = c(25, 80), breaks = seq(30, 70, 10)) +
  labs(
    title    = "The pocket does not enlarge, wherever the ligand sits",
    subtitle = paste("Every published substrate- or detergent-bound MexB",
                     "structure, measured in one common frame."),
    # plotmath rather than literal Å/³, so the labels render on any device
    # and in any locale
    x = expression("depth into the pocket from the periplasmic entrance ("
                   * ring(A) * ")"),
    y = expression("ligand-free pocket volume (" * ring(A)^3 * ")"),
    caption = paste0(
      "r = ", sprintf("%+.2f", r_depth), " between volume and depth; ",
      sprintf("%+.2f", r_size), " between volume and ligand size.\n",
      "Depth is arc length back from the periplasmic mouth along the widest ",
      "ligand-free entry channel of the reference protomer, to the ligand ",
      "centroid.\nEach protomer superposed on 39 pocket-lining C-alpha of one ",
      "reference, so the measuring sphere sits identically in every ",
      "structure. Engineered MexB chimeras excluded.")) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "#dde5e8"),
    axis.line          = element_line(colour = "#9fb0b8"),
    plot.title         = element_text(colour = "#104862", face = "bold",
                                      size = 14),
    plot.subtitle      = element_text(colour = "#4a5a66"),
    plot.caption       = element_text(colour = "#4a5a66", hjust = 0, size = 8),
    plot.caption.position = "plot",
    plot.title.position   = "plot")

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(aes(label = ligand), size = 3.4,
                                    seed = 1, min.segment.length = Inf,
                                    box.padding = 0.6, show.legend = FALSE)
} else {
  message("ggrepel not installed - using fixed label offsets")
  p <- p + geom_text(aes(label = ligand), size = 3.4, vjust = -1.3,
                     show.legend = FALSE)
}

out <- file.path(repo, "results", "figures", "poster")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out, "P4_R.pdf"), p, width = 9, height = 5.5)
ggsave(file.path(out, "P4_R.png"), p, width = 9, height = 5.5, dpi = 300)
cat("wrote", file.path(out, "P4_R.pdf"), "and .png\n")
