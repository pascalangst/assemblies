library(venn)
library(ggplot2)
library(ggpolypath)
library(RColorBrewer)
library(data.table)


pfam.reduced <- read.csv("~/Library/CloudStorage/OneDrive-UniversitätBasel/assemblies/data/pfam.reduced.csv")
pfam.reduced <- pfam.reduced[,c(1,3,5,6,4,8,7,2,9)]
x <- pfam.reduced[c(2:8)]
x[x> 0] <- 1
x <- x[!apply(x, 1, function(row) all(row == 0)), ]

colnames(x) <- gsub("\\.\\.", ". ", colnames(x))

venn(x, ilabels = "counts", zcolor = brewer.pal(7, "Set2"), col = brewer.pal(7, "Set2"), opacity = .75, box = FALSE, ilcs = 1, sncs = 1.5, lwd = 3)

y <- pfam.reduced

total_abundance <- rowSums(y[,c(2:8)])
y[order(total_abundance, decreasing = T),]

y <- y[apply(y[,c(2:8)], 1, function(row) any(row == 0)), ]
y$sum <- rowSums(y[,c(2:8)])
y[order(y$sum, decreasing = T),]

colnames(y) <- gsub("\\.\\.", ". ", colnames(y))
y_long <- melt(head(y[order(y$sum, decreasing = T),], n = 50), 
               id.vars = c("PFAM", "descriptions"), 
               measure.vars = c("M. daphniae", "G. daphniae", "C. obtusa", "B. daphniae", 
                                "H. tvaerminnensis", "O. colligata", "G. intestinalis", "sum"),
               variable.name = "species",
               value.name = "abundance")

y_long$PFAM <- factor(y_long$PFAM, rev(head(y[order(y$sum, decreasing = T),], n = 50)$PFAM))

ggplot(y_long, aes(x = PFAM, y = abundance, fill = species, width=0.75)) + 
  labs(x = "PFAM domain", y = "Abundance", fill = NULL) +
  geom_bar(stat = "identity") +
  facet_grid(cols = vars(species)) +
  coord_flip() +
  theme_bw() + 
  scale_fill_brewer(palette="Set2") +
  scale_x_discrete(labels=paste0(c(rev(y_long$descriptions)), " (", rev(y_long$PFAM), ")")) +
  theme(axis.text.x =element_text(size=15), axis.text.y =element_text(size=10), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_blank(),  
        legend.position = "none", strip.text = element_text(face = "italic"), strip.text.x = element_text(size = 10))
  
