# Swisspond pool of D. pulex infected with Gurleya daphniae (Graz)

# nextdenovo assembly (shasta is less complete and more fragmented)

# polishing with illumina and long reads 
bwa-mem2 index nd.asm.fasta
bwa-mem2 mem -t 8 -M nd.asm.fasta  BSSE_QGF_227929_HM7Y3DSX5_2_CH_H_Gurleya_CTGTGGCG_CACTAGCC_S3_L002_R1_001_MM_1.fastq.gz BSSE_QGF_227929_HM7Y3DSX5_2_CH_H_Gurleya_CTGTGGCG_CACTAGCC_S3_L002_R2_001_MM_1.fastq.gz > illumina.nd.asm.magna.mapped.sam
samtools view -Sb illumina.nd.asm.magna.mapped.sam | samtools sort -o - | tee illumina.nd.asm.mapped.bam | bedtools genomecov -ibam - > illumina.nd.asm.mapped.bam.txt
samtools depth -a illumina.nd.asm.mapped.bam |  awk '{sum+=$3} END { print "Average of  = ",sum/NR}' > illumina.nd.asm.mapped.bam.cov
samtools flagstat illumina.nd.asm.mapped.bam > illumina.nd.asm.mapped.bam.flagstat

conda activate nanopolish
ls CH-H-Gv.fastq > lgs.fofn
ls BSSE_QGF_227929_HM7Y3DSX5_2_CH_H_Gurleya_CTGTGGCG_CACTAGCC_S3_L002_R1_001_MM_1.fastq.gz BSSE_QGF_227929_HM7Y3DSX5_2_CH_H_Gurleya_CTGTGGCG_CACTAGCC_S3_L002_R2_001_MM_1.fastq.gz > sgs.fofn
nextPolish run_polish_sgs-lgs.cfg 

# polishing based on ont error profile
conda activate medaka
medaka_consensus -i CH-H-Gv.fastq -d polished_sgs-lgs/genome.nextpolish.fasta -o nextDenovo_polished_sgs-lgs_medaka -t 6 -m r941_min_high_g303

# polishing based on illumina reads
conda activate phd-Basics
cd nextDenovo_polished_sgs-lgs_medaka
bwa-mem2 index consensus.fasta
bwa-mem2 mem -t 10 consensus.fasta BSSE_QGF_227929_HM7Y3DSX5_2_CH_H_Gurleya_CTGTGGCG_CACTAGCC_S3_L002_R1_001_MM_1.fastq.gz BSSE_QGF_227929_HM7Y3DSX5_2_CH_H_Gurleya_CTGTGGCG_CACTAGCC_S3_L002_R2_001_MM_1.fastq.gz > illumina.consensus.sam
samtools view -@ 5 -Sb illumina.consensus.sam | samtools sort -@ 5 -T xyz -o illumina.consensus.bam
rm illumina.consensus.sam
samtools index -b illumina.consensus.bam

conda activate medaka
pilon --genome consensus.fasta --frags illumina.consensus.bam --changes --output consensus.pilon -Xmx110G

# check coverage per contig
bbmap.sh in=../BSSE_QGF_227929_HM7Y3DSX5_2_CH_H_Gurleya_CTGTGGCG_CACTAGCC_S3_L002_R1_001_MM_1.fastq.gz  in2=../BSSE_QGF_227929_HM7Y3DSX5_2_CH_H_Gurleya_CTGTGGCG_CACTAGCC_S3_L002_R2_001_MM_1.fastq.gz ref=consensus.pilon.fasta covstats=covstats.txt t=8

# purge haplotigs
conda activate phd-Basics
bwa-mem2 index consensus.pilon.fasta
bwa-mem2 mem -p -t 8 consensus.pilon.fasta ../BSSE_QGF_227929_HM7Y3DSX5_2_CH_H_Gurleya_CTGTGGCG_CACTAGCC_S3_L002_R1_001_MM_1.fastq.gz ../BSSE_QGF_227929_HM7Y3DSX5_2_CH_H_Gurleya_CTGTGGCG_CACTAGCC_S3_L002_R2_001_MM_1.fastq.gz > CH_H_Gurleya.sam
samtools view -Sb -@ 4 CH_H_Gurleya.sam | samtools sort -@ 4 -o CH_H_Gurleya.bam

bedtools genomecov -ibam CH_H_Gurleya.bam > CH_H_Gurleya.bam.txt

conda activate purge_haplotigs
purge_haplotigs  hist  -b CH_H_Gurleya.bam  -g consensus.pilon.fasta  -t 8
purge_haplotigs  cov  -i CH_H_Gurleya.bam.txt  -l 88  -m 176  -h 300  \
            [-o coverage_stats.csv -j 80  -s 80 ]
purge_haplotigs  purge  -g consensus.pilon.fasta  -c coverage_stats.csv -t 10

# check kmer distribution
minimap2 -ax map-ont -t 3 curated.fasta ../CH-H-Gv.fastq > curated.self.ont.sam
samtools flagstat curated.self.ont.sam > curated.self.ont.sam.flagstat.log
samtools view -Sb curated.self.ont.sam | samtools sort -o - | tee curated.self.ont.bam | bedtools genomecov -ibam - > curated.self.ont.bam.txt
samtools depth -a curated.self.ont.bam |  awk '{sum+=$3} END { print "Average of  = ",sum/NR}';

# make a fastq file out of a BAM
samtools view -b -F 4 curated.self.ont.bam > curated.self.ont.mapped.bam
bedtools bamtofastq -i curated.self.ont.mapped.bam -fq curated.self.ont.mapped.fq

# create input for the kmer histograms
conda activate medaka
jellyfish count -C -m 21 -s 1000000000 -t 10 curated.self.ont.mapped.fq -o curated.self.ont.mapped.jf
jellyfish histo -t 10 curated.self.ont.mapped.jf > curated.self.ont.mapped.histo
genomescope2 -i curated.self.ont.mapped.histo -o . -p 2 -k 21


# short reads:
bwa-mem2 index curated.fasta
bwa-mem2 mem -t 12 -M curated.fasta ../BSSE_QGF_227929_HM7Y3DSX5_2_CH_H_Gurleya_CTGTGGCG_CACTAGCC_S3_L002_R1_001_MM_1.fastq.gz ../BSSE_QGF_227929_HM7Y3DSX5_2_CH_H_Gurleya_CTGTGGCG_CACTAGCC_S3_L002_R2_001_MM_1.fastq.gz > curated.self.Illumina.sam
samtools flagstat curated.self.Illumina.sam > curated.self.Illumina.sam.flagstat.log
#273038270 + 0 in total (QC-passed reads + QC-failed reads)
#1853222 + 0 secondary
#0 + 0 supplementary
#0 + 0 duplicates
#61871163 + 0 mapped (22.66% : N/A)
#271185048 + 0 paired in sequencing
#135592524 + 0 read1
#135592524 + 0 read2
#57958164 + 0 properly paired (21.37% : N/A)
#58868594 + 0 with itself and mate mapped
#1149347 + 0 singletons (0.42% : N/A)
#705128 + 0 with mate mapped to a different chr
#413695 + 0 with mate mapped to a different chr (mapQ>=5)

samtools view -Sb curated.self.Illumina.sam | samtools sort -o - | tee curated.self.Illumina.bam | bedtools genomecov -ibam - > curated.self.Illumina.bam.txt
samtools depth -a curated.self.Illumina.bam |  awk '{sum+=$3} END { print "Average of  = ",sum/NR}';
#Average of = 448.226


# checking kmer distribution
# make a fastq file out of a BAM
samtools view -b -F 4 curated.self.Illumina.bam > curated.self.Illumina.mapped.bam
bedtools bamtofastq -i curated.self.Illumina.mapped.bam -fq curated.self.Illumina.mapped.fq

# create input for the kmer histograms
kmercountexact.sh in=curated.self.Illumina.mapped.fq khist=curated.self.Illumina.mapped.hist.txt peaks=curated.self.Illumina.mapped.peaks.txt t=12






conda activate basics
picard AddOrReplaceReadGroups \
    I=curated.self.Illumina.mapped.bam \
    O=curated.self.Illumina.mapped.rg.bam \
    RGPL=Illumina \
    RGLB=Lane1 \
    RGSM=BE-OM-3 \
    RGPU=TTAGGC

picard MarkDuplicates \
    INPUT=curated.self.Illumina.mapped.rg.bam OUTPUT=curated.self.Illumina.mapped.md.bam METRICS_FILE=curated.self.Illumina.mapped.rg.bam.txt -Xmx100G -XX:ParallelGCThreads=3
    
samtools index curated.self.Illumina.mapped.md.bam

picard CreateSequenceDictionary R=curated.fasta O=curated.dict
samtools faidx curated.fasta


gatk3 \
    -T HaplotypeCaller \
    -R curated.fasta \
    -I curated.self.Illumina.mapped.md.bam \
    -o curated.self.Illumina.g.vcf.gz \
    --emitRefConfidence GVCF
