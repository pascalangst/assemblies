library(venn)
library(ggplot2)
library(ggpolypath)
library(RColorBrewer)


pfam.reduced <- read.csv("~/Library/CloudStorage/OneDrive-UniversitätBasel/assemblies/data/pfam.reduced.csv")
x <- pfam.reduced[,c(3,5,6,4,8,7,2)]
x[x> 0] <- 1
x <- x[!apply(x, 1, function(row) all(row == 0)), ]

colnames(x) <- gsub("\\.\\.", ". ", colnames(x))

venn(x, ilabels = "counts", zcolor = brewer.pal(7, "Set1"), col = brewer.pal(7, "Set1"), opacity = .4, box = FALSE, ilcs = 1, sncs = 1.5)
