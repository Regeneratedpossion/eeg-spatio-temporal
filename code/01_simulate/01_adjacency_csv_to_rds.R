# Converts adjacency_61ch.csv (from make_adjacency_61ch.py) to the matrix stored in
# resources/adjacency/adjacency_xproj_yproj_61ch.rds.
m <- as.matrix(read.csv("adjacency_61ch.csv", row.names = 1, check.names = FALSE))
storage.mode(m) <- "integer"
saveRDS(m, "adjacency_xproj_yproj_61ch.rds")
