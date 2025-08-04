library(venn)
library(ggplot2)
library(RColorBrewer)
library(data.table)


pfam <- read.csv("~/Library/CloudStorage/OneDrive-UniversitätBasel/assemblies/data/funannotate_compare_Vnec_GiHiFi/pfam/pfam.results.csv")
pfam.reduced <- pfam[,c(1,2,3,7,25,26,27,28,29,30,31)]
x <- pfam.reduced[c(2:10)]
x[x> 0] <- 1
x <- x[!apply(x, 1, function(row) all(row == 0)), ]

colnames(x) <- gsub("\\.\\.", ". ", colnames(x))

venn(x, ilabels = "counts", zcolor = brewer.pal(7, "Set2"), col = brewer.pal(7, "Set2"), opacity = .75, box = FALSE, ilcs = 1, sncs = 1.5, lwd = 3)

y <- pfam.reduced

total_abundance <- rowSums(y[,c(2:10)])
y[order(total_abundance, decreasing = T),]

y <- y[apply(y[,c(2:10)], 1, function(row) any(row == 0)), ]
y$sum <- rowSums(y[,c(2:10)])
y$prop <- 1-rowSums(y[,c(2:10)] == 0)/ncol(y[,c(2:10)])
y[order(y$sum, decreasing = T),]

colnames(y) <- gsub("\\.\\.", ". ", colnames(y))
y_long <- melt(head(y[order(y$sum, decreasing = T),], n = 61), 
               id.vars = c("X", "descriptions"), 
               measure.vars = c("M. daphniae", "G. vavrai", "C. obtusa", "B. daphniae", 
                                "H. tvaerminnensis", "O. colligata", "E. intestinalis", "V. necatrix", "G. intestinalis", "sum", "prop"),
               variable.name = "species",
               value.name = "abundance")

y_long$X <- factor(y_long$X, rev(head(y[order(y$sum, decreasing = T),], n = 61)$X))
y_long$species <- sub("^(.{4}).*", "\\1.", y_long$species)
y_long$species <- factor(y_long$species, c("M. d.", "B. d.", "C. o.", "G. v.", "H. t.", "O. c.", "E. i.", "V. n.", "G. i.", "sum", "prop."))

ggplot(y_long, aes(x = X, y = abundance, fill = species, width=0.75)) + 
  labs(x = "PFAM domain", y = "Abundance", fill = NULL) +
  geom_bar(stat = "identity", fill = "black") +
  facet_grid(cols = vars(species)) +
  coord_flip() +
  theme_bw() + 
  #scale_fill_brewer(palette="Set3") +
  scale_x_discrete(labels=paste0(c(rev(y_long$descriptions)), " (", rev(y_long$X), ")")) +
  scale_y_continuous(breaks = c(0,40)) +
  theme(axis.text.x =element_text(size=15), axis.text.y =element_text(size=14), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_blank(),  
        legend.position = "none", strip.text = element_text(face = "italic"), strip.text.x = element_text(size = 10))

ggplot(y_long, aes(x = X, y = abundance, fill = species, width=0.75)) + 
  labs(x = "PFAM domain", y = "Abundance", fill = NULL) +
  geom_bar(stat = "identity", fill = "black") +
  facet_grid(cols = vars(species), scales = "free") +
  coord_flip() +
  theme_bw() + 
  #scale_fill_brewer(palette="Set3") +
  scale_x_discrete(labels=paste0(c(rev(y_long$descriptions)), " (", rev(y_long$X), ")")) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
  theme(axis.text.x =element_text(size=15), axis.text.y =element_text(size=14), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_blank(),  
        legend.position = "none", strip.text = element_text(face = "italic"), strip.text.x = element_text(size = 10))
