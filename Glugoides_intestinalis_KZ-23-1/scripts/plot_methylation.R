library(ggplot2)
library(ggpubr)
library(gsubfn)

modified_bases.5mC.bed <- read.delim("analysis/KZ-23-1.manual-merge.modified_bases.5mC.bed", header=FALSE)

my_files <- list.files(pattern = "tsv$", recursive=T, path = "analysis/nucleotide_composition_KZ-23-1.manual-merge/")
my_names <- gsub("KZ-23-1.manual.merge.", "", gsub(".tsv", "", my_files))
my_files <- paste0("analysis/nucleotide_composition_KZ-23-1.manual-merge/", my_files)
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
  
  p1 <- ggplot(na.omit(data), aes(V2, V11)) + geom_line() + ylim(0,100) + xlab(paste("")) + ylab("% 5-mC methyl. based on reads") + theme_bw()
  p3 <- ggplot(data_nucl, aes(V1, V2)) + geom_line() + ylim(min(data_nucl$V2),max(data_nucl$V2)) + xlab(paste("")) + ylab("% GC") + theme_bw()
  p4 <- ggplot(data, aes(V2, V5)) + geom_line() + ylim(0,100) + xlab(paste("position in", i)) + ylab("Number of reads") + theme_bw()
  arrange <- ggarrange(p1,p3,p4, ncol = 1, nrow = 3)
  print(arrange)
  ggsave(paste0(i, ".png"), arrange, height = 7.5, width = 7)
}




dict <- list("contig_01" = "scf1", 
             "contig_02" = "scf2", 
             "contig_03" = "scf3",
             "contig_04" = "scf4", 
             "contig_05" = "scf5", 
             "contig_06" = "scf6", 
             "contig_07" = "scf7", 
             "contig_08" = "scf9", 
             "contig_09" = "scf10", 
             "contig_10" = "scf11", 
             "contig_11" = "scf12", 
             "contig_12" = "scf13", 
             "contig_13" = "scf8") 

modified_bases.5mC.bed.meth <- modified_bases.5mC.bed[modified_bases.5mC.bed$V5 > 3 & modified_bases.5mC.bed$V11 > 50,]

homer.peak.file <- data.frame("Peak.ID"=rownames(modified_bases.5mC.bed.meth), "chromosome"=gsubfn("\\S+", dict, modified_bases.5mC.bed.meth$V1), "start"=modified_bases.5mC.bed.meth$V2, "end"=modified_bases.5mC.bed.meth$V3, "strand"=modified_bases.5mC.bed.meth$V6)

write.table(homer.peak.file, file="KZ-23-1.homer_50.peaks", sep = "\t", quote = F, row.names = F, col.names = F)

