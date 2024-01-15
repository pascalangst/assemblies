library(ggplot2)
library(ggpubr)

modified_bases.5mC.bed <- read.delim("~/Downloads/output.hifi.call_mods.modbam.pbmm2.freq.aggregate.all.bedmethyl", header=FALSE)

my_files <- list.files(pattern = "tsv$", recursive=T, path = "~/Downloads/nucleotide_composition_GB-ELK1-1.busco.blast.self/")
my_names <- gsub("GB-ELK1-1.busco.blast.self.", "", gsub(".tsv", "", my_files))
my_files <- paste0("~/Downloads/nucleotide_composition_GB-ELK1-1.busco.blast.self/", my_files)
my_tables <- lapply(my_files, read.table, header=F)
names(my_tables) <- my_names
for (i in 1:(length(my_tables))) {
  my_tables[[i]]$name <- names(my_tables[i])
}
nucleotide_composition <- as.data.frame(do.call("rbind",my_tables))

for (i in unique(modified_bases.5mC.bed$V1)){
  data <- modified_bases.5mC.bed[modified_bases.5mC.bed$V1==i,]
  data_nucl <- nucleotide_composition[nucleotide_composition$name==i,]
  
  p1 <- ggplot(data, aes(V2, V11)) + geom_point() + ylim(0,100) + xlab(paste("position in", i)) + ylab("% 5-mC methyl. based on reads") + theme_bw()
  #p2 <- ggplot(data, aes(V2, V12)) + geom_point() + ylim(0,100) + xlab(paste("position in", i)) + ylab("% 5-hmC methyl. based on reads") + theme_bw()
  p3 <- ggplot(data_nucl, aes(V1, V3)) + geom_point() + geom_smooth() + ylim(0,100) + xlab(paste("position in", i)) + ylab("% AT") + theme_bw()
  p4 <- ggplot(data, aes(V2, V5)) + geom_point() + ylim(0,80) + xlab(paste("position in", i)) + ylab("Number of reads") + geom_smooth() + theme_bw()
  #arrange <- ggarrange(p1,p2,p3,p4, ncol = 1, nrow = 4)
  arrange <- ggarrange(p1,p3,p4, ncol = 1, nrow = 3)
  #ggsave(paste0(i, ".png"), arrange, height = 10, width = 7)
  ggsave(paste0(i, ".png"), arrange, height = 7.5, width = 7)
}

