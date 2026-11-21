args <- commandArgs(trailingOnly = TRUE)
labels <- readRDS(args[[1]])
stopifnot(nrow(labels) == 341328L,
          identical(names(labels), c("Barcode", "Layer1", "Layer2", "Layer3")),
          !anyDuplicated(labels$Barcode))
primary <- labels[labels$Layer1 == "Primary Colon", ]
stopifnot(nrow(primary) == 60801L)
write.csv(primary, args[[2]], row.names = FALSE, fileEncoding = "UTF-8")
