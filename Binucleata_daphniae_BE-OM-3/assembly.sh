## BE-OM-3 infected with Binucleata daphniae

## tried many different approaches, Masurca was the best:

# trim adaptors
fastp -i /scicore/projects/openbis/userstore/duw_ebert/20230213185142181-60920231/BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R1_001_MM_1.fastq.gz -I /scicore/projects/openbis/userstore/duw_ebert/20230213185142181-60920231/BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R2_001_MM_1.fastq.gz -o BE_OM_3_Micro_2_2b.filtered.R1.fq.gz -O BE_OM_3_Micro_2_2b.filtered.R2.fq.gz  -q 30 -l 100 -w 55

# map to magna
bwa-mem2 mem -t 55 -M /scicore/home/ebertd/angpas00/PhD/data/ref/29082016_Xinb3_ref.fasta BE_OM_3_Micro_2_2b.filtered.R1.fq.gz BE_OM_3_Micro_2_2b.filtered.R2.fq.gz > BE_OM_3_Micro_2_2b.filtered.magna.mapped.sam

samtools view -@ 10 -Sb BE_OM_3_Micro_2_2b.filtered.magna.mapped.sam | samtools sort -@ 10 -o BE_OM_3_Micro_2_2b.filtered.magna.mapped.bam -

samtools view -@ 10 -Sb -f 4 BE_OM_3_Micro_2_2b.filtered.magna.mapped.bam | \
samtools sort -@ 10 -n -o BE_OM_3_Micro_2_2b.filtered.unmapped.sorted.bam -
samtools fastq BE_OM_3_Micro_2_2b.filtered.unmapped.sorted.bam -1 Bd.filtered_reads_R1.fastq.gz -2 Bd.filtered_reads_R2.fastq.gz -0 /dev/null -s /dev/null -n

# masurca assembly
masurca -t 16 -i Bd.filtered_reads_R1.fastq.gz,Bd.filtered_reads_R2.fastq.gz -r no_magna/Bd_reads.fastq -p masurca

cd CA.mr.49.17.15.0.02/

# check for contamination
~/PhD/data/trimmed_reads/spring_2017/SK-44_spr2017.megahit_asm/spades/Misc/runTaxonomizedBLAST.pl -t 12 -query primary.genome.scf.fasta
bbmap.sh in=../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R1_001_MM_1.fastq.gz in2=../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R2_001_MM_1.fastq.gz ref=primary.genome.scf.fasta covstats=covstats.masurca.txt t=12

# purge haplotigs
bwa-mem2 index primary.genome.scf.fasta
bwa-mem2 mem -p -t 8 primary.genome.scf.fasta ../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R1_001_MM_1.fastq.gz ../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R2_001_MM_1.fastq.gz > primary.genome.scf.sam
samtools view -Sb -@ 4 primary.genome.scf.sam | samtools sort -@ 4 -o primary.genome.scf.bam
bedtools genomecov -ibam primary.genome.scf.bam > primary.genome.scf.bam.txt

conda activate purge_haplotigs
purge_haplotigs  hist  -b primary.genome.scf.bam  -g primary.genome.scf.fasta  -t 8
purge_haplotigs  cov  -i primary.genome.scf.bam.txt  -l 125  -m 220  -h 400              [-o coverage_stats.csv -j 80  -s 80 ]
purge_haplotigs  purge  -g primary.genome.scf.fasta  -c coverage_stats.csv -t 10

# remove contaminants
seqtk subseq curated.fasta <(grep -v "scf7180000000611\|scf7180000000666\|scf7180000000681\|scf7180000000684" <(cut -f1 covstats.masurca.txt)) > curated.nomagna.fasta

# purge haplotigs again
bwa-mem2 index curated.nomagna.fasta
bwa-mem2 mem -p -t 16 curated.nomagna.fasta ../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R1_001_MM_1.fastq.gz ../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R2_001_MM_1.fastq.gz > curated.nomagna.sam
samtools view -Sb -@ 4 curated.nomagna.sam | samtools sort -@ 4 -o curated.nomagna.bam
bedtools genomecov -ibam curated.nomagna.bam > curated.nomagna.bam.txt

purge_haplotigs  hist  -b curated.nomagna.bam  -g curated.nomagna.fasta  -t 8
purge_haplotigs  cov  -i curated.nomagna.bam.txt  -l 125  -m 220  -h 400              [-o coverage_stats.csv -j 80  -s 80 ]
purge_haplotigs  purge  -g curated.nomagna.fasta  -c coverage_stats.csv -t 10 -o curated2
