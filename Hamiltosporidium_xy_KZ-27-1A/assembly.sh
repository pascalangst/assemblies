## KZ-27 (dead) infected with some Hamiltosporidium spp. (somewhere between tvaer and magni)
# second sample (spores) is almost exclusively contaminants

### mapping to magna
bwa-mem2 mem -t 8 -M ~/PhD/data/ref/29082016_Xinb3_ref.fasta BSSE_QGF_228027_H7LVWDSX5_3_KZ_27_1A_CTCATCAC_GGTGGCAC_S57_L003_R1_001_MM_1.fastq.gz BSSE_QGF_228027_H7LVWDSX5_3_KZ_27_1A_CTCATCAC_GGTGGCAC_S57_L003_R2_001_MM_1.fastq.gz > KZ_27.magna.mapped.sam
bwa-mem2 mem -t 8 -M ~/PhD/data/ref/29082016_Xinb3_ref.fasta BSSE_QGF_228028_H7LVWDSX5_3_KZ_27_1A_spores_ACCGATTA_GCCGTGGC_S58_L003_R1_001_MM_1.fastq.gz BSSE_QGF_228028_H7LVWDSX5_3_KZ_27_1A_spores_ACCGATTA_GCCGTGGC_S58_L003_R2_001_MM_1.fastq.gz > KZ-27_spores.magna.mapped.sam

samtools view -Sb -f 4 KZ_27.magna.mapped.sam | \
samtools sort -n -o KZ_27.unmapped.sorted.bam -
samtools fastq KZ_27.unmapped.sorted.bam -1 UP_reads_R1.fastq.gz -2 UP_reads_R2.fastq.gz -0 /dev/null -s /dev/null -n

samtools view -Sb KZ_27.magna.mapped.sam | samtools sort -o - | tee KZ_27.magna.mapped.bam | bedtools genomecov -ibam - > KZ_27.magna.mapped.bam.txt
samtools depth -a KZ_27.magna.mapped.bam |  awk '{sum+=$3} END { print "Average of  = ",sum/NR}';


# shovill assembly of reads not mapping to host
shovill --outdir shovill --R1 UP_reads_R1.fastq.gz --R2 UP_reads_R2.fastq.gz --cpus 12 --ram 64 --trim


# megahit assembly of reads not mapping to host (better!)
megahit -1 UP_reads_R1.fastq.gz  -2 UP_reads_R2.fastq.gz  -t 8  -o megahit

# get coverage per contig: filters are cov > 50x, length >= 500 bp, and GC < 0.33
bbmap.sh in=../UP_reads_R1.fastq.gz in2=../UP_reads_R2.fastq.gz ref=final.contigs.fa covstats=covstats.txt t=55
awk '{if ($5 > 50 && $6 > 499 && $7 < 0.33) print $1}' covstats.txt > contigs.cov50.500.GC33.names 
seqtk subseq final.contigs.fa contigs.cov50.500.GC33.names > contigs.cov50.500.GC33.fasta

# remove potential pulex contigs
bbmap.sh in=../FI-SK-58-2c.R1.fq in2=../FI-SK-58-2c.R2.fq  ref=contigs.cov50.500.GC33.fasta covstats=covstats.cured.txt t=55
awk '{if ($5 < 1) print $0}' covstats.cured.txt > nocured.names
seqtk subseq contigs.cov50.500.GC33.fasta nocured.names > contigs.cov50.500.GC33.nocured.fasta

# check for other contaminants
~/PhD/data/trimmed_reads/spring_2017/SK-44_spr2017.megahit_asm/spades/Misc/runTaxonomizedBLAST.pl -t 8 -query contigs.cov50.500.GC33.nocured.fasta
