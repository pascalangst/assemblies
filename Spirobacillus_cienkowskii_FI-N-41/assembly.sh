## FI-N-41 infected with Spirobacillus

# mapping to magna
bwa-mem2 mem -t 55 -M /scicore/home/ebertd/angpas00/PhD/data/ref/29082016_Xinb3_ref.fasta /scicore/projects/openbis/userstore/duw_ebert/20230213182836462-60920221/BSSE_QGF_227932_HM7Y3DSX5_2_FI_N_41_Sc_GCCGCAAC_TCATGTCT_S6_L002_R1_001_MM_1.fastq.gz /scicore/projects/openbis/userstore/duw_ebert/20230213182836462-60920221/BSSE_QGF_227932_HM7Y3DSX5_2_FI_N_41_Sc_GCCGCAAC_TCATGTCT_S6_L002_R2_001_MM_1.fastq.gz > FI_N_41_Sc.magna.mapped.sam

samtools view -Sb FI_N_41_Sc.magna.mapped.sam | samtools sort -o - | tee FI_N_41_Sc.magna.mapped.bam | bedtools genomecov -ibam - > FI_N_41_Sc.magna.mapped.bam.txt
samtools depth -a FI_N_41_Sc.magna.mapped.bam |  awk '{sum+=$3} END { print "Average of  = ",sum/NR}' > FI_N_41_Sc.magna.mapped.bam.cov

samtools view -Sb -f 4 FI_N_41_Sc.magna.mapped.bam | \
samtools sort -n -o FI_N_41_Sc.unmapped.sorted.bam -
samtools fastq FI_N_41_Sc.unmapped.sorted.bam -1 Sc_reads_R1.fastq.gz -2 Sc_reads_R2.fastq.gz -0 /dev/null -s /dev/null -n

# shovill assembly of reads not mapping to host
shovill --outdir shovill --R1 Sc_reads_R1.fastq.gz --R2 Sc_reads_R2.fastq.gz --cpus 12 --ram 64


# unicycler assembly of reads not mapping to host (better!)
unicycler -1 Sc2.filtered_reads_R1.0.1.fastq.gz -2 Sc2.filtered_reads_R2.0.1.fastq.gz -o unicycler_downsample_0.1 -t 55
cd unicycler_downsample_0.1

# get coverage per contig for infected and non-infected samples
bbmap.sh in=FI-SK-58-2c.R1.fq in2=FI-SK-58-2c.R2.fq  ref=assembly.fasta covstats=covstats.cured.txt t=55
bbmap.sh in=../Sc2.filtered_reads_R1.fastq.gz in2=Sc2.filtered_reads_R2.fastq.gz ref=assembly.fasta covstats=covstats.txt t=55

# Two versions: 
#FI_N_41_Sc.FULL.contigs.fasta	FULL version
#FI_N_41_Sc.contigs.fasta		reduced based on sequence length and possible contaminant contigs





