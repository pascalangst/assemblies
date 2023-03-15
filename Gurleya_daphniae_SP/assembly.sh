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
