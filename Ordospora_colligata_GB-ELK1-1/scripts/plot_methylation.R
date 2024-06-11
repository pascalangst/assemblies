library(ggplot2)
library(ggpubr)

modified_bases.5mC.bed <- read.delim("analysis/output.hifi.call_mods.modbam.pbmm2.freq.aggregate.all.manual_merge.bed", header=FALSE)

my_files <- list.files(pattern = "tsv$", recursive=T, path = "analysis/nucleotide_composition_GB-LK1-1.manual.merge/")
my_names <- gsub("GB-LK1-1.manual.merge.", "", gsub(".tsv", "", my_files))
my_files <- paste0("analysis/nucleotide_composition_GB-LK1-1.manual.merge/", my_files)
my_tables <- lapply(my_files, read.table, header=F)
names(my_tables) <- my_names
for (i in 1:(length(my_tables))) {
  my_tables[[i]]$name <- names(my_tables[i])
}
nucleotide_composition <- as.data.frame(do.call("rbind",my_tables))

for (i in unique(modified_bases.5mC.bed$V1)){
  data <- modified_bases.5mC.bed[modified_bases.5mC.bed$V1==i,]
  data <- data[data$V5 > 5,]
  data_nucl <- nucleotide_composition[nucleotide_composition$name==i,]
  
  data <- data %>%
    mutate(V11 = zoo::rollapply(data$V11, width = 20, by = 10, FUN = mean, fill = NA, align = "right"))
  
  p1 <- ggplot(na.omit(data), aes(V2, V11)) + geom_line() + ylim(0,10) + xlab(paste("")) + ylab("% 5-mC methyl. based on reads") + theme_bw()
  p3 <- ggplot(data_nucl, aes(V1, V2)) + geom_line() + ylim(min(data_nucl$V2),max(data_nucl$V2)) + xlab(paste("")) + ylab("% GC") + theme_bw()
  p4 <- ggplot(data, aes(V2, V5)) + geom_line() + ylim(0,100) + xlab(paste("position in", i)) + ylab("Number of reads") + theme_bw()
  arrange <- ggarrange(p1,p3,p4, ncol = 1, nrow = 3)
  print(arrange)
  ggsave(paste0(i, ".png"), arrange, height = 10, width = 7)
}

