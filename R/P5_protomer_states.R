# Pocket volume in every protomer of every MexB structure, by conformational
# state. Self-contained: paste into R and run. Needs ggplot2.
#
# 39 protomers from 9 structures (four carry two trimers in the asymmetric
# unit). Free volume in one 16 A sphere fixed in the frame of the reference
# Binding protomer, with every protomer superposed on its 39 pocket-lining
# C-alpha, so the sphere sits at the same anatomical position throughout.
# In an Access or Extrusion protomer that is not the protomer's own pocket as
# it would be defined in isolation - it is how open the substrate site is at
# the same place.
#
# Full table, including the four sphere radii and the connected-volume
# measure: results/tables/protomer_pockets.csv

library(ggplot2)

d <- read.csv(text = "pdb,chain,state,volume,ours
Amp_MexB_20260826,D,Binding,1561,TRUE
Amp_MexB_20260826,E,Binding,1969,TRUE
Amp_MexB_20260826,F,Extrusion,1638,TRUE
MexB_DDM_3_20260730,D,Access,1158,TRUE
MexB_DDM_3_20260730,E,Binding,2102,TRUE
MexB_DDM_3_20260730,F,Extrusion,1734,TRUE
21FO,A,Access,1102,FALSE
21FO,B,Binding,2148,FALSE
21FO,C,Extrusion,1722,FALSE
21FP,A,Access,1138,FALSE
21FP,B,Binding,2233,FALSE
21FP,C,Extrusion,1804,FALSE
2V50,A,Access,1425,FALSE
2V50,B,Binding,2036,FALSE
2V50,C,Extrusion,1604,FALSE
2V50,D,Access,1265,FALSE
2V50,E,Binding,2741,FALSE
2V50,F,Extrusion,1604,FALSE
3W9I,A,Access,1068,FALSE
3W9I,B,Binding,2350,FALSE
3W9I,C,Extrusion,1620,FALSE
3W9I,D,Access,1055,FALSE
3W9I,E,Binding,2146,FALSE
3W9I,F,Extrusion,1738,FALSE
3W9J,A,Access,1052,FALSE
3W9J,B,Binding,1973,FALSE
3W9J,C,Extrusion,1760,FALSE
3W9J,D,Access,1018,FALSE
3W9J,E,Binding,1890,FALSE
3W9J,F,Extrusion,1738,FALSE
6IIA,A,Access,1044,FALSE
6IIA,B,Binding,2109,FALSE
6IIA,C,Extrusion,1612,FALSE
6IIA,D,Access,1046,FALSE
6IIA,E,Binding,2045,FALSE
6IIA,F,Extrusion,1719,FALSE
6T7S,J,Binding,1607,FALSE
6T7S,K,Access,1342,FALSE
6T7S,L,Access,1512,FALSE
", stringsAsFactors = FALSE)

d$state <- factor(d$state, levels = c("Access", "Binding", "Extrusion"))

# the poster's protomer colours, darkened to stay legible on white
pal <- c(Access = "#0A9DA0", Binding = "#CA0FC1", Extrusion = "#0F9C1B")

aggregate(volume ~ state, d, function(v)
  c(n = length(v), mean = mean(v), sd = sd(v)))

set.seed(0)
ggplot(d, aes(state, volume, colour = state)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.62,
               colour = "#16202a", linewidth = 0.5) +
  geom_jitter(data = subset(d, !ours), width = 0.16, height = 0,
              size = 3, alpha = 0.55) +
  geom_jitter(data = subset(d, ours), width = 0.16, height = 0,
              size = 4.5, shape = 18) +
  scale_colour_manual(values = pal, guide = "none") +
  labs(x = NULL,
       y = expression("free volume at the substrate site (" * ring(A)^3 * ")")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.line          = element_line(colour = "#9fb0b8"))
