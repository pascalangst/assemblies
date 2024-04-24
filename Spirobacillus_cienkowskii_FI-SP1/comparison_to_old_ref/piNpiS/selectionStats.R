library(ggplot2)

selectionStats_Spiro <- read.delim("~/OneDrive - Universität Basel/assemblies/Spirobacillus_cienkowskii_FI-SP1/comparison_to_old_ref/piNpiS/selectionStats_Spiro.txt")
selectionStats_Spiro$Alignment <- gsub(".aligned.fas.oneline", "", selectionStats_Spiro$Alignment)

selectionStats_Spiro$pin <- selectionStats_Spiro$PiN / selectionStats_Spiro$NumN
selectionStats_Spiro$pis <- selectionStats_Spiro$PiS / selectionStats_Spiro$NumS

selectionStats_Spiro$pi_ratio <- selectionStats_Spiro$pin /selectionStats_Spiro$pis

summary(selectionStats_Spiro)
hist(selectionStats_Spiro$pin, breaks = 20)
hist(selectionStats_Spiro$pis, breaks = 20)
hist(selectionStats_Spiro$pi_ratio, breaks = 20)

buscos <- read.delim("~/OneDrive - Universität Basel/assemblies/Spirobacillus_cienkowskii_FI-SP1/comparison_to_old_ref/piNpiS/full_table.tsv", header=FALSE, comment.char="#")

selectionStats_Spiro$type <- NA
selectionStats_Spiro[selectionStats_Spiro$Alignment %in% buscos$V3, ]$type <- "BUSCO"
selectionStats_Spiro[selectionStats_Spiro$Alignment %in% c("Spiro2_04520", "Spiro2_04605", "Spiro2_02780"), ]$type <- "Chitinase"
selectionStats_Spiro[selectionStats_Spiro$Alignment %in% c("Spiro2_02390", "Spiro2_05315", "Spiro2_09515","Spiro2_09515", "Spiro2_12215","Spiro2_03345","Spiro2_03330","Spiro2_03340","Spiro2_03335","Spiro2_03350","Spiro2_08300","Spiro2_08300","Spiro2_08300","Spiro2_01195","Spiro2_01200"), ]$type <- "Carotenoid pathway" 
selectionStats_Spiro[selectionStats_Spiro$Alignment %in% c("Spiro2_09460", "Spiro2_05065", "Spiro2_02510", "Spiro2_08925", "Spiro2_05235","Spiro2_08105", "Spiro2_03885","Spiro2_10955", "Spiro2_04390","Spiro2_08255", "Spiro2_04255", "Spiro2_04260"), ]$type <- "Glycolysis"

# BUSCOs in other groups are only in the other group, not in BUSCO
ggplot(selectionStats_Spiro, aes(x=type, y=pi_ratio)) + geom_boxplot() + theme_bw() +
  theme(legend.position="top", legend.box="vertical", axis.text=element_text(size=15), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), legend.text = element_text(size=15)) +
  labs(x=paste0(""),
       y=bquote(pi[N]/pi[S]))

ggplot(selectionStats_Spiro, aes(x=type, y=Pi)) + geom_boxplot() + theme_bw() +
  theme(legend.position="top", legend.box="vertical", axis.text=element_text(size=15), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), legend.text = element_text(size=15)) +
  labs(x=paste0(""),
       y=bquote(pi~axis~cut)) + ylim(c(0,0.1))
