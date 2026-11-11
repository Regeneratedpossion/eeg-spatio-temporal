# Thesis Figure 2.3 (fig:electrode_adjacency).
source(here::here("code/00_config.R")); setwd(thesis_dir)   # reads and writes relative to the thesis folder
# Draws the 61-electrode adjacency graph used for spatial adjacency.
suppressPackageStartupMessages(library(ggplot2))
x <- readLines("fig_data/adjacency_edges.txt")
k <- which(x == "##")
el <- read.csv(text = x[1:(k - 1)])
ed <- read.csv(text = x[(k + 1):length(x)])
ed <- merge(merge(ed, el, by.x = "from", by.y = "channel"), el, by.x = "to", by.y = "channel",
            suffixes = c("1", "2"))
added <- data.frame(from = c("AF1","AF2","AF7","AF8","P3","PO1","P4","PO2","O1","O2"),
                    to = c("AF7","AF8","F3","F4","PO7","PO7","PO8","PO8","POZ","POZ"))
key <- function(a, b) paste(pmin(a, b), pmax(a, b))
ed$added <- key(ed$from, ed$to) %in% key(added$from, added$to)
stopifnot(sum(ed$added) == 10, nrow(ed) == 200)
p <- ggplot() +
  geom_segment(data = ed[!ed$added, ], aes(x = xproj1, y = yproj1, xend = xproj2, yend = yproj2), colour = "grey55", linewidth = 0.4) +
  geom_segment(data = ed[ed$added, ], aes(x = xproj1, y = yproj1, xend = xproj2, yend = yproj2), colour = "#C0392B", linewidth = 0.6, linetype = "dashed") +
  geom_point(data = el, aes(xproj, yproj), size = 5.5, colour = "grey30", fill = "white", shape = 21) +
  geom_text(data = el, aes(xproj, yproj, label = channel), size = 1.9) +
  coord_equal() + theme_void()
ggsave(file.path(fig_dir(thesis_path("diagnostics")), "electrode_adjacency.pdf"), p, width = 5.6, height = 5.6)
