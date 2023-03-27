library(ggplot2)
library(ggpubr)

#meth.path <- "~/Library/CloudStorage/OneDrive-UniversitätBasel/assemblies/Binucleata_daphniae_BE-OM-3/methylation/Binucleata_daphniae_modified_bases.5mC.bed"
#nucl.path <- "~/Library/CloudStorage/OneDrive-UniversitätBasel/assemblies/Binucleata_daphniae_BE-OM-3/nucl_bias/curated2/"

meth.path <- "~/Library/CloudStorage/OneDrive-UniversitätBasel/assemblies/Gurleya_daphniae_SP/methylation/Gurleya_daphniae_modified_bases.5mC.bed"
nucl.path <- "~/Library/CloudStorage/OneDrive-UniversitätBasel/assemblies/Gurleya_daphniae_SP/nucl_bias/curated/"

# Methylation information from Megalodon
modified_bases.5mC.bed <- read.delim(meth.path, header=FALSE)

# Nucleotide composition
my_files <- list.files(pattern = "tsv$", recursive=T, path = nucl.path)
my_names <- gsub(paste0(basename(nucl.path),"."), "", gsub(".tsv", "", my_files))
my_files <- paste0(nucl.path, my_files)
my_tables <- lapply(my_files, read.table, header=F)
names(my_tables) <- my_names
for (i in 1:(length(my_tables))) {
  my_tables[[i]]$name <- names(my_tables[i])
}
nucleotide_composition <- as.data.frame(do.call("rbind",my_tables))

# plot methylation and read coverage
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


