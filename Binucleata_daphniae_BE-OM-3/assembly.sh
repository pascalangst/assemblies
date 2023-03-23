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



# checking kmer distribution

# long reads:
minimap2 -ax map-ont -t 8 curated2.fasta BE-OM-3_depl_pass.fastq > curated2.self.ont.sam
samtools flagstat curated2.self.ont.sam > curated2.self.ont.sam.flagstat.log
samtools view -Sb curated2.self.ont.sam | samtools sort -o - | tee curated2.self.ont.bam | bedtools genomecov -ibam - > curated2.self.ont.bam.txt
samtools depth -a curated2.self.ont.bam |  awk '{sum+=$3} END { print "Average of  = ",sum/NR}' > curated2.self.ont.bam.cov
samtools view -b -F 4 curated2.self.ont.bam > curated2.self.ont.mapped.bam
bedtools bamtofastq -i curated2.self.ont.mapped.bam -fq curated2.self.ont.mapped.fq

# create input for the kmer histograms
conda activate medaka
jellyfish count -C -m 21 -s 1000000000 -t 10 curated2.self.ont.mapped.fq -o curated2.self.ont.mapped.jf
jellyfish histo -t 10 curated2.self.ont.mapped.jf > curated2.self.ont.mapped.histo
genomescope2 -i curated2.self.ont.mapped.histo -o . -p 2 -k 21 -o diplo
genomescope2 -i curated2.self.ont.mapped.histo -o . -p 4 -k 21 -o tetra

# short reads:
bwa-mem2 index curated2.fasta
bwa-mem2 mem -t 12 -M curated2.fasta ../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R1_001_MM_1.fastq.gz ../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R2_001_MM_1.fastq.gz > curated2.self.Illumina.sam
samtools flagstat curated2.self.Illumina.sam > curated2.self.Illumina.sam.flagstat.log
#190339675 + 0 in total (QC-passed reads + QC-failed reads)
#345971 + 0 secondary
#0 + 0 supplementary
#0 + 0 duplicates
#43191025 + 0 mapped (22.69% : N/A)
#189993704 + 0 paired in sequencing
#94996852 + 0 read1
#94996852 + 0 read2
#40780028 + 0 properly paired (21.46% : N/A)
#42182262 + 0 with itself and mate mapped
#662792 + 0 singletons (0.35% : N/A)
#1232314 + 0 with mate mapped to a different chr
#592052 + 0 with mate mapped to a different chr (mapQ>=5)

samtools view -Sb curated2.self.Illumina.sam | samtools sort -o - | tee curated2.self.Illumina.bam | bedtools genomecov -ibam - > curated2.self.Illumina.bam.txt
samtools depth -a curated2.self.Illumina.bam |  awk '{sum+=$3} END { print "Average of  = ",sum/NR}';
#Average of = 619.332


# checking kmer distribution
# make a fastq file out of a BAM
samtools view -b -F 4 curated2.self.Illumina.bam > curated2.self.Illumina.mapped.bam
bedtools bamtofastq -i curated2.self.Illumina.mapped.bam -fq curated2.self.Illumina.mapped.fq

# create input for the kmer histograms
kmercountexact.sh in=curated2.self.Illumina.mapped.fq khist=curated2.self.Illumina.mapped.hist.txt peaks=curated2.self.Illumina.mapped.peaks.txt t=12







picard AddOrReplaceReadGroups \
    I=curated2.self.Illumina.mapped.bam \
    O=curated2.self.Illumina.mapped.rg.bam \
    RGPL=Illumina \
    RGLB=Lane1 \
    RGSM=BE-OM-3 \
    RGPU=TTAGGC

picard MarkDuplicates \
    INPUT=curated2.self.Illumina.mapped.rg.bam OUTPUT=curated2.self.Illumina.mapped.md.bam METRICS_FILE=curated2.self.Illumina.mapped.rg.bam.txt -Xmx100G -XX:ParallelGCThreads=3
    
samtools index curated2.self.Illumina.mapped.md.bam

picard CreateSequenceDictionary R=curated2.fasta O=curated2.dict
samtools faidx curated2.fasta


gatk3 \
    -T HaplotypeCaller \
    -R curated2.fasta \
    -I curated2.self.Illumina.mapped.md.bam \
    -o curated2.self.Illumina.g.vcf.gz \
    --emitRefConfidence GVCF
    
    
    



## nextdenovo assembly
# polishing with illumina and long reads
conda activate nanopolish
nextPolish run_polish_sgs-lgs.cfg

# polishing based on ont error profile
conda activate mamba-medaka
medaka_consensus -i BE-OM-3_depl_pass.fastq -d polished_sgs-lgs/genome.nextpolish.fasta -o nextDenovo_polished_sgs-lgs_medaka -t 6 -m r941_min_high_g303

# polishing based on illumina reads
conda activate phd-Basics
cd nextDenovo_polished_sgs-lgs_medaka
bwa-mem2 index consensus.fasta
bwa-mem2 mem -t 10 consensus.fasta ../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R1_001_MM_1.fastq.gz ../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R2_001_MM_1.fastq.gz > illumina.consensus.sam
samtools view -@ 5 -Sb illumina.consensus.sam | samtools sort -@ 5 -T xyz -o illumina.consensus.bam
rm illumina.consensus.sam
samtools index -b illumina.consensus.bam

conda activate mamba-medaka
pilon --genome consensus.fasta --frags illumina.consensus.bam --changes --output consensus.pilon -Xmx110G



## shasta assembly
# polishing with illumina and long reads
conda activate nanopolish
nextPolish run_polish_sgs-lgs.cfg
# polishing based on ont error profile
conda activate mamba-medaka
medaka_consensus -i ../BE-OM-3_depl_pass.fastq -d polished_sgs-lgs/genome.nextpolish.fasta -o shasta_polished_sgs-lgs_medaka -t 6 -m r941_min_high_g303

# polishing based on illumina reads
conda activate phd-Basics
cd shasta_polished_sgs-lgs_medaka
bwa-mem2 index consensus.fasta
bwa-mem2 mem -t 10 consensus.fasta ../../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R1_001_MM_1.fastq.gz ../../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R2_001_MM_1.fastq.gz > illumina.consensus.sam
samtools view -@ 5 -Sb illumina.consensus.sam | samtools sort -@ 5 -T xyz -o illumina.consensus.bam
rm illumina.consensus.sam
samtools index -b illumina.consensus.bam

conda activate mamba-medaka
pilon --genome consensus.fasta --frags illumina.consensus.bam --changes --output consensus.pilon -Xmx110G




## merge shasta and nextdenovo
merge_wrapper.py nextDenovo_polished_sgs-lgs_medaka/consensus.pilon.fasta shasta/shasta_polished_sgs-lgs_medaka/consensus.pilon.fasta -pre nd_shasta_pilon
merge_wrapper.py shasta/shasta_polished_sgs-lgs_medaka/consensus.pilon.fasta nextDenovo_polished_sgs-lgs_medaka/consensus.pilon.fasta -pre shasta_nd_pilon
# merge more complete with more contiguous assembly
merge_wrapper.py shasta/shasta_polished_sgs-lgs_medaka/consensus.pilon.fasta merged_nd_shasta_pilon.fasta -pre shasta_merge_nd_shasta_pilon



# purge haplotigs from nd assembly then merge
conda activate phd-Basics
cd nextDenovo_polished_sgs-lgs_medaka/
bwa-mem2 index consensus.pilon.fasta
bwa-mem2 mem -p -t 8 consensus.pilon.fasta ../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R1_001_MM_1.fastq.gz ../BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R2_001_MM_1.fastq.gz > nd.sam
samtools view -Sb -@ 4 nd.sam | samtools sort -@ 4 -o nd.bam

bedtools genomecov -ibam nd.bam > nd.bam.txt

conda activate purge_haplotigs
purge_haplotigs  hist  -b nd.bam  -g consensus.pilon.fasta  -t 8
purge_haplotigs  cov  -i nd.bam.txt  -l 125  -m 225  -h 400  \
            [-o coverage_stats.csv -j 80  -s 80 ]
purge_haplotigs  purge  -g consensus.pilon.fasta  -c coverage_stats.csv -t 10

cd ../
conda activate quickmerge
merge_wrapper.py nextDenovo_polished_sgs-lgs_medaka/curated.fasta shasta/shasta_polished_sgs-lgs_medaka/consensus.pilon.fasta -pre nd-curated_shasta_pilon
#merge_wrapper.py shasta/shasta_polished_sgs-lgs_medaka/consensus.pilon.fasta nextDenovo_polished_sgs-lgs_medaka/consensus.pilon.fasta -pre shasta_nd_pilon
# merge more complete with more contiguous assembly
merge_wrapper.py shasta/shasta_polished_sgs-lgs_medaka/consensus.pilon.fasta merged_nd-curated_shasta_pilon.fasta -pre shasta_merge_nd-curated_shasta_pilon

# purge haplotigs from merged assembly
conda activate phd-Basics

bwa-mem2 index merged_shasta_merge_nd-curated_shasta_pilon.fasta
bwa-mem2 mem -p -t 8 merged_shasta_merge_nd-curated_shasta_pilon.fasta BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R1_001_MM_1.fastq.gz BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R2_001_MM_1.fastq.gz > merged_shasta_merge_nd-curated_shasta_pilon.sam
samtools view -Sb -@ 4 merged_shasta_merge_nd-curated_shasta_pilon.sam | samtools sort -@ 4 -o merged_shasta_merge_nd-curated_shasta_pilon.bam

bedtools genomecov -ibam merged_shasta_merge_nd-curated_shasta_pilon.bam > merged_shasta_merge_nd-curated_shasta_pilon.bam.txt

conda activate purge_haplotigs
purge_haplotigs  hist  -b merged_shasta_merge_nd-curated_shasta_pilon.bam  -g merged_shasta_merge_nd-curated_shasta_pilon.fasta  -t 8
purge_haplotigs  cov  -i merged_shasta_merge_nd-curated_shasta_pilon.bam.txt  -l 125  -m 225  -h 400  \
            [-o coverage_stats.csv -j 80  -s 80 ]
purge_haplotigs  purge  -g merged_shasta_merge_nd-curated_shasta_pilon.fasta  -c coverage_stats.csv -t 10

bbmap.sh in=BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R1_001_MM_1.fastq.gz  in2=BSSE_QGF_227942_HM7Y3DSX5_2_BE_OM_3_Micro_2_2b_AGCATGGA_AGCGCTAA_S16_L002_R2_001_MM_1.fastq.gz ref=curated.fasta covstats=covstats.curated.txt t=16
~/PhD/data/trimmed_reads/spring_2017/SK-44_spr2017.megahit_asm/spades/Misc/runTaxonomizedBLAST.pl -t 8 -query curated.fasta

seqtk subseq curated.fasta <(grep -v "770_np1212_pilon\|724_np1212_pilon\|894_np1212_pilon\|734_np1212_pilon\|474_np1212_pilon" <(cut -f1 covstats.curated.txt)) > curated.nomagna.fasta


minimap2 -ax map-ont -t 3 curated.nomagna.fasta BE-OM-3_depl_pass.fastq > curated.nomagna.self.ont.sam
samtools flagstat curated.nomagna.self.ont.sam > curated.nomagna.self.ont.sam.flagstat.log
samtools view -Sb curated.nomagna.self.ont.sam | samtools sort -o - | tee curated.nomagna.self.ont.bam | bedtools genomecov -ibam - > curated.nomagna.self.ont.bam.txt
samtools depth -a curated.nomagna.self.ont.bam |  awk '{sum+=$3} END { print "Average of  = ",sum/NR}' > curated.nomagna.self.ont.bam.cov

# checking kmer distribution
# make a fastq file out of a BAM
samtools view -b -F 4 curated.nomagna.self.ont.bam > curated.nomagna.self.ont.mapped.bam
bedtools bamtofastq -i curated.nomagna.self.ont.mapped.bam -fq curated.nomagna.self.ont.mapped.fq

# create input for the kmer histograms
conda activate medaka
jellyfish count -C -m 21 -s 1000000000 -t 10 curated.nomagna.self.ont.mapped.fq -o curated.nomagna.self.ont.mapped.jf
jellyfish histo -t 10 curated.nomagna.self.ont.mapped.jf > curated.nomagna.self.ont.mapped.histo
genomescope2 -i curated.nomagna.self.ont.mapped.histo -o . -p 2 -k 21 -o diplo
genomescope2 -i curated.nomagna.self.ont.mapped.histo -o . -p 4 -k 21 -o tetra
