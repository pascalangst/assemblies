# GB-ELK1-1 infected with Ordospora colligata (same isolate as used for previous reference assembly)
module load minimap2/2.20-GCCcore-10.3.0
minimap2 -ax map-hifi -t 12 29082016_Xinb3_ref.fasta GBELK11_Circular_Consensus_Sequencing_Reads/m64156_230929_124402.bc2002--bc2002.hifi_reads.fastq.gz > GB-ELK1-1_Xinb3.sam

module load SAMtools/1.7-goolf-1.7.20
samtools view -Sb -f 4 GB-ELK1-1_Xinb3.sam | \
samtools sort -n -o GB-ELK1-1_Xinb3.hifi.unmapped.bam -
samtools fastq GB-ELK1-1_Xinb3.hifi.unmapped.bam -0 GB-ELK1-1_Xinb3.hifi.unmapped.fastq.gz -n

# assembly using hifiasm (also tried fly assembly and meta assembly with the same read set)
hifiasm -o no-host/GB-ELK1-1.no-host.asm -t8 GB-ELK1-1_Xinb3.hifi.unmapped.fastq.gz

# identify Oc contigs

# get contigs with microsporidia BUSCOs
busco -c 12 -m geno -i GB-ELK1-1.no-host.asm.bp.p_ctg.fa -o BUSCO5.FULL.output -l microsporidia_odb10
awk '{print $3}' BUSCO5.FULL.output/run_microsporidia_odb10/full_table.tsv | sort | uniq > contigs.busco.names
seqtk subseq GB-ELK1-1.no-host.asm.bp.p_ctg.fa contigs.busco.names > GB-ELK1-1.busco.fa

# remove remaining contaminants
~/PhD/data/trimmed_reads/spring_2017/SK-44_spr2017.megahit_asm/spades/Misc/runTaxonomizedBLAST.pl -t 12 -query GB-ELK1-1.busco.fa
seqtk subseq GB-ELK1-1.no-host.asm.bp.p_ctg.fa contigs.busco.names > GB-ELK1-1.busco.blast.fa
seqtk subseq GB-ELK1-1.no-host.asm.bp.p_ctg.fa contigs.busco.names > GB-ELK1-1.busco.blast.self.fa

# identify rRNA genes
makeblastdb -in GB-ELK1-1.busco.blast.self.fa -dbtype nucl
blastn -query ~/master_thesis/ordospora/reference/16S.fasta  -db GB-ELK1-1.busco.blast.self.fa -max_target_seqs 100 -out blast.16S.out -outfmt 6

# investigate regions with HGT
minimap2 -ax map-hifi -t 12 GB-ELK1-1.busco.blast.self.fa GB-ELK1-1_Xinb3.hifi.unmapped.fastq.gz > GB-ELK1-1_GB-ELK1-1.nohost.sam
minimap2 -ax map-hifi -t 12 GB-ELK1-1.busco.blast.self.fa m64156_230929_124402.bc2002--bc2002.hifi_reads.fastq.gz > GB-ELK1-1_GB-ELK1-1.sam
samtools view -@6  -Sb GB-ELK1-1_GB-ELK1-1.nohost.sam | samtools sort -@6 -o GB-ELK1-1_GB-ELK1-1.nohost.bam 
samtools depth GB-ELK1-1_GB-ELK1-1.nohost.bam > GB-ELK1-1_GB-ELK1-1.nohost.depth

# plot in R
#GB.ELK1.1_GB.ELK1.1.nohost <- read.delim("~/Downloads/GB-ELK1-1_GB-ELK1-1.nohost.depth", header=FALSE)
#plot(x=GB.ELK1.1_GB.ELK1.1.nohost[GB.ELK1.1_GB.ELK1.1.nohost$V1 == "ptg000012l",]$V2, y=GB.ELK1.1_GB.ELK1.1.nohost[GB.ELK1.1_GB.ELK1.1.nohost$V1 == "ptg000012l",]$V3)
#plot(x=GB.ELK1.1_GB.ELK1.1.nohost[GB.ELK1.1_GB.ELK1.1.nohost$V1 == "ptg000012l",]$V2, y=GB.ELK1.1_GB.ELK1.1.nohost[GB.ELK1.1_GB.ELK1.1.nohost$V1 == "ptg000012l",]$V3, type = "l")
#plot(x=GB.ELK1.1_GB.ELK1.1.nohost[GB.ELK1.1_GB.ELK1.1.nohost$V1 == "ptg000003l",]$V2, y=GB.ELK1.1_GB.ELK1.1.nohost[GB.ELK1.1_GB.ELK1.1.nohost$V1 == "ptg000003l",]$V3, type = "l")
#plot(x=GB.ELK1.1_GB.ELK1.1.nohost[GB.ELK1.1_GB.ELK1.1.nohost$V1 == "ptg000015l",]$V2, y=GB.ELK1.1_GB.ELK1.1.nohost[GB.ELK1.1_GB.ELK1.1.nohost$V1 == "ptg000015l",]$V3, type = "l")

#GB.ELK1.1_GB.ELK1.1 <- read.delim("~/Downloads/GB-ELK1-1_GB-ELK1-1.depth", header=FALSE)
#plot(x=GB.ELK1.1_GB.ELK1.1[GB.ELK1.1_GB.ELK1.1$V1 == "ptg000012l",]$V2, y=GB.ELK1.1_GB.ELK1.1[GB.ELK1.1_GB.ELK1.1$V1 == "ptg000012l",]$V3, type = "l")
#plot(x=GB.ELK1.1_GB.ELK1.1.nohost[GB.ELK1.1_GB.ELK1.1.nohost$V1 == "ptg000012l",]$V2, y=GB.ELK1.1_GB.ELK1.1.nohost[GB.ELK1.1_GB.ELK1.1.nohost$V1 == "ptg000012l",]$V3, type = "l")
