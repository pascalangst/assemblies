## SK-39 infected with Larssonia

### mapping to magna
##bwa-mem2 mem -t 8 -M ~/PhD/data/ref/29082016_Xinb3_ref.fasta BSSE_QGF_228026_H7LVWDSX5_3_SK_39_ACACATTC_CTGTCGGT_S56_L003_R1_001_MM_1.fastq.gz BSSE_QGF_228026_H7LVWDSX5_3_SK_39_ACACATTC_CTGTCGGT_S56_L003_R2_001_MM_1.fastq.gz > SK-39_Larssonia.magna.mapped.sam
##
##samtools view -Sb -f 4 SK-39_Larssonia.magna.mapped.sam | \
##samtools sort -n -o SK-39_Larssonia.unmapped.sorted.bam -
##samtools fastq SK-39_Larssonia.unmapped.sorted.bam -1 Lo_reads_R1.fastq.gz -2 Lo_reads_R2.fastq.gz -0 /dev/null -s /dev/null -n
##
##shovill --outdir shovill --R1 Lo_reads_R1.fastq.gz --R2 Lo_reads_R2.fastq.gz --cpus 12 --ram 64 --trim
##
##samtools view -Sb SK-39_Larssonia.magna.mapped.sam | samtools sort -o - | tee SK-39_Larssonia.magna.mapped.bam | bedtools genomecov -ibam - > SK-39_Larssonia.magna.mapped.bam.txt
##samtools depth -a SK-39_Larssonia.magna.mapped.bam |  awk '{sum+=$3} END { print "Average of  = ",sum/NR}';

# it is not magna




### mapping to longispina
##ln -s  /home/ebertgrp/assemblies/daphnia_longispina/pacbio/hifi/D_longispina_inb4.v0.1.fasta.gz . 
##
##bwa-mem2 index D_longispina_inb4.v0.1.fasta.gz
##
##bwa-mem2 mem -t 8 -M D_longispina_inb4.v0.1.fasta.gz BSSE_QGF_228026_H7LVWDSX5_3_SK_39_ACACATTC_CTGTCGGT_S56_L003_R1_001_MM_1.fastq.gz BSSE_QGF_228026_H7LVWDSX5_3_SK_39_ACACATTC_CTGTCGGT_S56_L003_R2_001_MM_1.fastq.gz > SK-39_Larssonia.longispina.mapped.sam
##
##samtools view -Sb SK-39_Larssonia.longispina.mapped.sam | samtools sort -o - | tee SK-39_Larssonia.longispina.mapped.bam | bedtools genomecov -ibam - > SK-39_Larssonia.longispina.mapped.bam.txt
##samtools depth -a SK-39_Larssonia.longispina.mapped.bam |  awk '{sum+=$3} END { print "Average of  = ",sum/NR}' > SK-39_Larssonia.longispina.mapped.bam.cov
##
##samtools flagstat SK-39_Larssonia.longispina.mapped.bam > SK-39_Larssonia.longispina.mapped.bam.flagstat

# it is not longispina




### mapping to pulex
bwa-mem2 index GCA_911175335.1_PA42_4.2_genomic.fna.gz

bwa-mem2 mem -t 8 -M GCA_911175335.1_PA42_4.2_genomic.fna.gz BSSE_QGF_228026_H7LVWDSX5_3_SK_39_ACACATTC_CTGTCGGT_S56_L003_R1_001_MM_1.fastq.gz BSSE_QGF_228026_H7LVWDSX5_3_SK_39_ACACATTC_CTGTCGGT_S56_L003_R2_001_MM_1.fastq.gz > SK-39_Larssonia.pulex.mapped.sam

samtools view -Sb SK-39_Larssonia.pulex.mapped.sam | samtools sort -o - | tee SK-39_Larssonia.pulex.mapped.bam | bedtools genomecov -ibam - > SK-39_Larssonia.pulex.mapped.bam.txt
samtools depth -a SK-39_Larssonia.pulex.mapped.bam |  awk '{sum+=$3} END { print "Average of  = ",sum/NR}';

samtools view -Sb -f 4 SK-39_Larssonia.pulex.mapped.bam | \
samtools sort -n -o SK-39_Larssonia.unmapped_pulex.sorted.bam -
samtools fastq SK-39_Larssonia.unmapped_pulex.sorted.bam -1 Lo_from_pulex_reads_R1.fastq.gz -2 Lo_from_pulex_reads_R2.fastq.gz -0 /dev/null -s /dev/null -n


# shovill assembly of reads not mapping to host
conda activate shovill

shovill --outdir shovill_from_pulex --R1 Lo_from_pulex_reads_R1.fastq.gz --R2 Lo_from_pulex_reads_R2.fastq.gz --cpus 12 --ram 64 --trim

# remove potential pulex contigs
cd shovill_from_pulex
bbmap.sh in=/home/pascal/PhD/data/trimmed_reads/spring_2017/SK-44_spr2017.megahit_asm/no_host/SRR10159717.fastq.gz ref=contigs.fa covstats=covstats.pulex.txt interleaved=T
awk '{if ($8 < 1) print $1}' covstats.pulex.txt  > contigs.pulex.names
seqtk subseq contigs.fa contigs.pulex.names > contigs.non-pulex.fa

# repeat with another sample (from close-by) 
bbmap.sh in=/home/pascal/PhD/data/trimmed_reads/summer_2017/SK-44_smr2017.pe.fq.gz ref=contigs.non-pulex.fa covstats=covstats.smr17.txt interleaved=T
awk '{if ($8 < 1) print $1}' covstats.smr17.txt  > contigs.smr17.names
seqtk subseq contigs.non-pulex.fa contigs.smr17.names > contigs.nosmr17.fa

# check for other contaminants
 ~/PhD/data/trimmed_reads/spring_2017/SK-44_spr2017.megahit_asm/spades/Misc/runTaxonomizedBLAST.pl -t 8 -query contigs.non-pulex.fa

# get coverage per contig
bbmap.sh in=../Lo_from_pulex_reads_R1.fastq.gz in2=../Lo_from_pulex_reads_R2.fastq.gz ref=contigs.nosmr17.fa covstats=covstats.nosmr17.txt t=55

# filters are cov > 50x, length >= 500 bp, and GC < 0.33
awk '{if ($8 > 50 && $9 > 499 && $10 < 0.33) print $1}' covstats.nosmr17.txt > contigs.cov50.500.GC33.noSMR.names 
seqtk subseq contigs.nosmr17.fa contigs.cov50.500.GC33.noSMR.names  > contigs.cov50.500.GC33.noSMR.fasta



# megahit assembly of reads not mapping to host (better!)
megahit -1 Lo_from_pulex_reads_R1.fastq.gz  -2 Lo_from_pulex_reads_R2.fastq.gz  -t 55  -o megahit
cd megahit

# get coverage per contig: filters are cov > 50x, length >= 500 bp, and GC < 0.33
bbmap.sh in=../Lo_from_pulex_reads_R1.fastq.gz in2=../Lo_from_pulex_reads_R2.fastq.gz ref=final.contigs.fa covstats=covstats.txt t=55
awk '{if ($5 > 50 && $6 > 499 && $7 < 0.33) print $1}' covstats.txt > contigs.cov50.500.GC33.names 
seqtk subseq final.contigs.fa contigs.cov50.500.GC33.names > contigs.cov50.500.GC33.fasta

# remove potential pulex contigs
bbmap.sh in=../SRR10159717.fastq.gz ref=contigs.cov50.500.GC33.fasta covstats=covstats.pulex.txt interleaved=T # no pulex detected

bbmap.sh in=../SK-44_smr2017.pe.fq.gz ref=contigs.cov50.500.GC33.fasta covstats=covstats.smr17.txt interleaved=T t=55
awk '{if ($5 < 1) print $0}' covstats.smr17.txt > no-smr17.names
seqtk subseq contigs.cov50.500.GC33.fasta no-smr17.names > contigs.cov50.500.GC33.noSMR.fasta

# check for other contaminants
~/PhD/data/trimmed_reads/spring_2017/SK-44_spr2017.megahit_asm/spades/Misc/runTaxonomizedBLAST.pl -t 8 -query contigs.cov50.500.GC33.noSMR.fasta

# checking kmer distribution:
bwa-mem2 index contigs.cov50.500.GC33.noSMR.fasta
bwa-mem2 mem -t 12 -M contigs.cov50.500.GC33.noSMR.fasta ../BSSE_QGF_228026_H7LVWDSX5_3_SK_39_ACACATTC_CTGTCGGT_S56_L003_R1_001_MM_1.fastq.gz ../BSSE_QGF_228026_H7LVWDSX5_3_SK_39_ACACATTC_CTGTCGGT_S56_L003_R2_001_MM_1.fastq.gz > contigs.cov50.500.GC33.noSMR.self.Illumina.sam
samtools flagstat contigs.cov50.500.GC33.noSMR.self.Illumina.sam > contigs.cov50.500.GC33.noSMR.self.Illumina.sam.flagstat.log


samtools view -Sb contigs.cov50.500.GC33.noSMR.self.Illumina.sam | samtools sort -o - | tee contigs.cov50.500.GC33.noSMR.self.Illumina.bam | bedtools genomecov -ibam - > contigs.cov50.500.GC33.noSMR.self.Illumina.bam.txt
samtools depth -a contigs.cov50.500.GC33.noSMR.self.Illumina.bam |  awk '{sum+=$3} END { print "Average of  = ",sum/NR}';


# checking kmer distribution
# make a fastq file out of a BAM
samtools view -b -F 4 contigs.cov50.500.GC33.noSMR.self.Illumina.bam > contigs.cov50.500.GC33.noSMR.self.Illumina.mapped.bam
bedtools bamtofastq -i contigs.cov50.500.GC33.noSMR.self.Illumina.mapped.bam -fq contigs.cov50.500.GC33.noSMR.self.Illumina.mapped.fq

# create input for the kmer histograms
kmercountexact.sh in=contigs.cov50.500.GC33.noSMR.self.Illumina.mapped.fq khist=contigs.cov50.500.GC33.noSMR.self.Illumina.mapped.hist.txt peaks=contigs.cov50.500.GC33.noSMR.self.Illumina.mapped.peaks.txt t=12


picard AddOrReplaceReadGroups \
    I=contigs.cov50.500.GC33.noSMR.self.Illumina.mapped.bam \
    O=contigs.cov50.500.GC33.noSMR.self.Illumina.mapped.rg.bam \
    RGPL=Illumina \
    RGLB=Lane1 \
    RGSM=Larssonia \
    RGPU=TTAGGC

picard MarkDuplicates \
    INPUT=contigs.cov50.500.GC33.noSMR.self.Illumina.mapped.rg.bam OUTPUT=contigs.cov50.500.GC33.noSMR.self.Illumina.mapped.md.bam METRICS_FILE=contigs.cov50.500.GC33.noSMR.self.Illumina.mapped.rg.bam.txt -Xmx100G -XX:ParallelGCThreads=3
    
samtools index contigs.cov50.500.GC33.noSMR.self.Illumina.mapped.md.bam

picard CreateSequenceDictionary R=contigs.cov50.500.GC33.noSMR.fasta O=contigs.cov50.500.GC33.noSMR.dict
samtools faidx contigs.cov50.500.GC33.noSMR.fasta


gatk3 \
    -T HaplotypeCaller \
    -R contigs.cov50.500.GC33.noSMR.fasta \
    -I contigs.cov50.500.GC33.noSMR.self.Illumina.mapped.md.bam \
    -o contigs.cov50.500.GC33.noSMR.self.Illumina.g.vcf.gz \
    --emitRefConfidence GVCF
