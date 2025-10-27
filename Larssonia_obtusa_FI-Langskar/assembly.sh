# base-calling (GPU) R10.4.1 flowcell V14 chemistry

/scicore/home/ebertd/angpas00/software/dorado/dorado-1.1.1-linux-x64/bin/dorado duplex /scicore/home/ebertd/angpas00/software/dorado/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 -v pod5/ --emit-fastq > duplex_251015.fastq

# align to longispina host reference and get unmapped reads
minimap2 -ax lr:hq -t 12 D_longispina_inb4.v0.1.fasta.gz duplex_251015.fastq.gz > Larssonia.inb4.ont.mapped.sam
samtools flagstat Larssonia.inb4.ont.mapped.sam > Larssonia.inb4.ont.mapped.sam.flagstat # 16.5 % unmapped
samtools view -Sb -f 4 Larssonia.inb4.ont.mapped.sam | \
samtools sort -n -o Larssonia.inb4.ont.unmapped.only.bam -
samtools fastq Larssonia.inb4.ont.unmapped.only.bam -0 Larssonia.inb4.ont.unmapped.only.fastq.gz -n

#flye --nano-hq duplex_251015.fastq.gz --meta --threads 55 --out-dir flye_ont_meta
#flye --nano-hq Larssonia.inb4.ont.unmapped.only.fastq.gz -g 18.13m --threads 55 --out-dir flye_ont_Larssonia
#flye --nano-hq Larssonia.inb4.ont.unmapped.only.fastq.gz --meta --threads 55 --out-dir flye_ont_Larssonia_meta
#~/bioinformatics/shasta-Linux-0.14.0 --config config --assemblyDirectory Shasta_1K_unmapped --input Larssonia.inb4.ont.unmapped.only.fastq 

# NextDenovo assembly (just unmapped reads) > 1KB
~/bioinformatics/NextDenovo/nextDenovo run.cfg 

## Assembly checks
# BUSCO score
busco -m geno -l ~/busco_downloads/lineages/microsporidia_odb12 -c 12 -i flye_ont_Larssonia/assembly.fasta -o BUSCO.flye_Larssonia.all
awk '{print $3}' BUSCO.Nextdenovo_1K.all/run_microsporidia_odb12/full_table.tsv | sort | uniq # get contigs with BUSCOS

# GC per contig
seqkit fx2tab -n -g 01_rundir/03.ctg_graph/nd.asm.fasta

# fcs (foreign contamination screen)
python3 ~/software/ncbi-fcs/fcs.py screen genome --fasta nd.asm.fasta --out-dir ./gx_out/ --gx-db ~/software/ncbi-fcs/gxdb --tax-id 2483382 
cat nd.asm.fasta | python3 ~/software/ncbi-fcs/fcs.py clean genome --action-report ./gx_out/nd.asm.2483382.fcs_gx_report.txt --output Larssonia.clean.fasta --contam-fasta-out nd.asm.contam.fasta

# ctg000380: GC outside others, and loop 
# ctg000090, ctg000120, ctg000290, ctg000300, ctg000380: no buscos and no fcs microsporidia hit

# polish (optional)
medaka_consensus -i Larssonia.inb4.ont.unmapped.only.fastq  -d 01_rundir/03.ctg_graph/nd.asm.fasta -o nextDenovo_medaka -t 18 -m r1041_e82_400bps_sup_v5.2.0

# purge haplotigs
cd nextDenovo_medaka
minimap2 -ax lr:hq -t 12 consensus.fasta ../Larssonia.inb4.ont.unmapped.only.fastq > Larssonia.consensus.ont.mapped.sam
samtools view -Sb -@ 4 Larssonia.consensus.ont.mapped.sam | samtools sort -@ 4 -o Larssonia.consensus.ont.mapped.bam
rm Larssonia.consensus.ont.mapped.sam
bedtools genomecov -ibam Larssonia.consensus.ont.mapped.bam > Larssonia.consensus.ont.mapped.bam.txt
conda activate purge_haplotigs
purge_haplotigs  hist  -b Larssonia.consensus.ont.mapped.bam  -g consensus.fasta  -t 8
purge_haplotigs  cov  -i Larssonia.consensus.ont.mapped.bam.txt  -l 50  -m 120  -h 160  \
            [-o coverage_stats.csv -j 80  -s 80 ]
purge_haplotigs  purge  -g consensus.fasta  -c coverage_stats.csv -t 10

# removed after purge_haplotigs
# ctg000380	-	-	-	-	JUNK
# ctg000120	-	-	-	-	JUNK
# ctg000290	-	-	-	-	JUNK
# ctg000090	-	-	-	-	JUNK
# 
# -> kept ctg000300 from above list


# NextDenovo assembly (all reads) > 1KB
~/bioinformatics/NextDenovo/nextDenovo run_all.cfg

# polish (optional)
medaka_consensus -i duplex_251015.fastq.gz -d 01_rundir_all/03.ctg_graph/nd.asm.fasta -o nextDenovo_medaka_all -t 6 -m r1041_e82_400bps_sup_v5.2.0


# merge more complete (nextdenovo with unmapped reads, medaka corrected, and purged haplotigs) with more contiguous assembly (nextdenovo from all reads)
merge_wrapper.py nextDenovo_medaka/curated.fasta nextDenovo_medaka_all/consensus.fasta -pre nextDenovo_medaka_nextdenovo_medaka_all

# check rRNAs
~/bioinformatics/infernal-1.1.2/src/cmscan -Z 36.350914 --cpu 2 --cut_ga --default --nohmmonly --fmt 2 --tblout rfam.tblout --clanin ~/bioinformatics/infernal-1.1.2/Rfam.clanin ~/bioinformatics/infernal-1.1.2/Rfam.cm merged_nextDenovo_medaka_nextdenovo_medaka_all.fasta  > rfam.cmscan
less rfam.tblout | grep Euk | awk '{print $4,$10,$11,$5,$5,$12}' | tr ' ' '\t' | awk '{OFS="\t"; if ($6=="+") {print} else {print $1,$3,$2,$4,$5,$6}}' | sort -k1,1 -k2,2n | bedtools merge -d 1000 -i - | awk '{OFS="\t"; {print "rRNA_"NR,$1,$2,$3,"+","rRNA"}}' 
 
# check telomeres
perl ../scripts/check_for_telomeres.pl -f merged_nextDenovo_medaka_nextdenovo_medaka_all.fasta -o TELOMERES.Nextdenovo_1K.medaka.curated.merge
