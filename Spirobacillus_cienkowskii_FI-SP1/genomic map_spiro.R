library(ComplexHeatmap)
library(stringr)
library(dplyr)
library(circlize)
library(data.table)
library(zoo)
library(tidyverse)

setwd("~/OneDrive - Universität Basel/assemblies/Spirobacillus_cienkowskii_FI-SP1/bakta_tmp/circos/")
#################################################################################
#### Extract CDS rRNA tRNA etc ####
################################################################################

spiro_prot_plus= read.delim("features-plus.txt", sep=" ", header=F,stringsAsFactors = FALSE) ### features plus strand 
spiro_prot_minus= read.delim("features-minus.txt", sep=" ", header=F,stringsAsFactors = FALSE)### features minus strand 
colnames(spiro_prot_minus)=c("contig_names","xstart","xend","orientation","color")
colnames(spiro_prot_plus)=c("contig_names","xstart","xend","orientation","color")
spiro_prot_plus=spiro_prot_plus %>% mutate( ystart= 2.1, yend=3) ### create the coordinates of the features  
spiro_prot_minus=spiro_prot_minus %>% mutate( ystart= 1, yend=1.9) ### create the coordinates of the features

spiro_size=2806830 

### there is probably an other way to read the color 
spiro_prot_plus=mutate(spiro_prot_plus ,col= case_when(color == "color=178,223,138" ~ "#B2DF8A", 
                                                       color == "color=204,204,204" ~ "#cccccc",
                                                       color == "color=251,128,114" ~ "#FB8072",
                                                       color == "color=253,180,98" ~ "#FDB462"))
spiro_prot_minus=mutate(spiro_prot_minus ,col= case_when(color == "color=178,223,138" ~ "#B2DF8A", 
                                                       color == "color=204,204,204" ~ "#cccccc",
                                                       color == "color=251,128,114" ~ "#FB8072",
                                                       color == "color=190,186,218" ~ "#BEBADA",
                                                       color == "color=253,180,98" ~ "#FDB462"))


### GC content ###
spiro_GC_content= read.delim("gc_content.txt", sep=" ", header=F,stringsAsFactors = FALSE)
colnames(spiro_GC_content)=c("contig_names","xstart","xend","value","color")
### you need to have one column for the pos value and one for the neg value to plot with 2 different colors
spiro_GC_content=mutate(spiro_GC_content ,value_pos= case_when(value >= 0 ~ value,value < 0  ~ 0))
spiro_GC_content=mutate(spiro_GC_content ,value_neg= case_when(value >= 0 ~ 0,value < 0  ~ value)) 
spiro_GC_content=spiro_GC_content[,c("contig_names","xstart","xend","value_pos","value_neg")]

### GC skew ###
spiro_GC_skew=read.delim("gc_skew.txt", sep=" ", header=F,stringsAsFactors = FALSE)
colnames(spiro_GC_skew)=c("contig_names","xstart","xend","value","color")
spiro_GC_skew=mutate(spiro_GC_skew ,value_pos= case_when(value >= 0 ~ value,value < 0  ~ 0)) 
spiro_GC_skew=mutate(spiro_GC_skew ,value_neg= case_when(value >= 0 ~ 0,value < 0  ~ value)) 
spiro_GC_skew=spiro_GC_skew[,c("contig_names","xstart","xend","value_pos","value_neg")]


### Plot ###
fs=15 ### font size 
### stupid way to create the axis legend
a=c(0)
b=0
for (i in 1:28){
  b=b+0.1
  a=c(a,b)
} 

a
brk <- a*10^6
### if you want to artificially increase the size of the gene: here doesn't work your genome is too dense 
add=400




pdf("try_6.pdf",width=10,height=10)
#tiff("try_5.tiff",,width=15,height=15,units="cm", res=800)
circos.clear()
df = data.frame(
  name  = c("contig_1"),
  start = c(0),
  end   = c(spiro_size))
circos.par(gap.degree=0, start.degree= 90 , cell.padding = c(0.01, 0, 0.01, 0))

circos.initialize(sectors=c("contig_1"),xlim=c(0,spiro_size))

# Genomic Features
# + strand
circos.track(ylim = c(1, 3), 
             
             #bg.col = c("#e0e0e0"), 
             bg.border = NA, track.height = 0.175, panel.fun = function(x, y){circos.axis(h="top",major.at=brk,labels=round(brk/10^6,1),labels.cex=0.9,lwd=0.7,labels.facing="clockwise")})
#circos.lines(c(0,1767095),c(1,1))
for (i in 1:nrow(spiro_prot_plus)){
  circos.rect(spiro_prot_plus$xstart[i]-add, spiro_prot_plus$ystart[i],
              spiro_prot_plus$xend[i]+add, spiro_prot_plus$yend[i], 
              col = spiro_prot_plus$col[i], border = spiro_prot_plus$col[i],lwd=0.01)
}
for (i in 1:nrow(spiro_prot_minus)){
  circos.rect(spiro_prot_minus$xstart[i]-add, spiro_prot_minus$ystart[i],
              spiro_prot_minus$xend[i]+add, spiro_prot_minus$yend[i], 
              col =  spiro_prot_minus$col[i], border =  spiro_prot_minus$col[i],lwd=0.01)
}
circos.text(0, 2, "A", cex=0.9,font = 2)
#GC content
circos.genomicTrack(data=spiro_GC_content,ylim=c(-0.205,0.20),bg.border = NA, track.height = 0.15,panel.fun=function(region,value,...) {
  circos.genomicLines(region,value,type="l",border=c("#33A02C","#E31A1C"),col=c("#33A02C","#E31A1C")
                      , lwd=0.3,baseline=c(0,0), area=TRUE, )
  circos.text(0, 0.12, "B", cex=0.9,font = 2)
})

# GC skew 
circos.genomicTrack(data=spiro_GC_skew,ylim=c(-0.26,0.35),bg.border = NA, track.height = 0.15,panel.fun=function(region,value,...) {
  circos.genomicLines(region,value,type="l",border=c("#FDBF6F","#1F78B4"),col=c("#FDBF6F","#1F78B4"), lwd=0.3, baseline=c(0,0), area=TRUE )
  circos.text(0, 0.35, "C", cex=0.9,font = 2)
})

lgd_points_A_1= Legend(labels=c("CDS", "rRNA"),labels_gp=gpar(fontsize = fs),
                       legend_gp =gpar(fill = c("#cccccc","#fb8072")), title_position = "lefttop", 
                       title = "A",title_gp = gpar(fontsize = fs, fontface = "bold"),grid_height = unit(4, "mm"), grid_width = unit(4.5, "mm"))
lgd_points_A_2= Legend(labels=c("tRNA/tmRNA", "ncRNA"),labels_gp=gpar(fontsize = fs),
                       legend_gp =gpar(fill = c("#b2df8a","#fdb462")), title_position = "lefttop", 
                       grid_height = unit(4, "mm"), grid_width = unit(4.5, "mm"))

lgd_points_C = Legend(labels=c("GC content above average", "GC content below average"),labels_gp=gpar(fontsize = fs), title_position = "lefttop",  #type = "lines",
                      legend_gp =gpar(fill=c("#33A02C","#E31A1C")),
                      title = "B",title_gp = gpar(fontsize =fs, fontface = "bold"),grid_height = unit(4, "mm"), grid_width = unit(4.5, "mm"))
lgd_points_D = Legend(labels=c("+ GC skew", "\u2212 GC skew"),labels_gp=gpar(fontsize = fs),title_position = "lefttop", #type = "lines",
                      legend_gp =gpar(fill=c("#FDBF6F","#1F78B4")),
                      title = "C",title_gp = gpar(fontsize = fs, fontface = "bold"),grid_height = unit(4, "mm"), grid_width = unit(4.5, "mm"))

pd_1 = packLegend(lgd_points_A_1,lgd_points_A_2 , direction="horizontal", gap = unit(1, "mm"))
pd_2 = packLegend(lgd_points_C, lgd_points_D, direction="vertical", gap = unit(1, "mm"))
pd=packLegend(pd_1, pd_2, direction = "vertical", gap = unit(1, "mm"))
draw(pd, just = c("left", "bottom"),x = unit(10, "cm"), y = unit(11, "cm"))


dev.off()

