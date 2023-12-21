# assembled using hifiasm v.0.16.0 without further parameter optimization

hifiasm -o longi.asm -t8 ../Pool1_Circular_Consensus_Sequencing_Reads/m64156_230621_203451.bc2096--bc2096.hifi_reads.fastq.gz
hifiasm -o denti.asm -t8 ../Pool1_Circular_Consensus_Sequencing_Reads/m64156_230621_203451.bc2095--bc2095.hifi_reads.fastq.gz


# to identify non-target contigs (same for both species):

# 1) Blast to latest nt database 

blastn -db /scicore/data/managed/BLAST_FASTA/latest/nt \
       -query longi.p_ctg.fa \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -out blast.scicore.out \
       -num_threads 64
       
# 2) Blast to latest uniprot database 
       
mkdir ../uniprot

cd ../uniprot

tar -xvf /scicore/data/managed/UniProt/latest/knowledgebase/reference_proteomes/Reference_Proteomes_2023_05.tar

touch reference_proteomes.fasta.gz
find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz

printf "accession\taccession.version\ttaxid\tgi\n" > reference_proteomes.taxid_map
zcat */*/*.idmapping.gz | grep "NCBI_TaxID" | awk '{print $1 "\t" $1 "\t" $3 "\t" 0}' >> reference_proteomes.taxid_map

diamond makedb -p 16 --in reference_proteomes.fasta.gz --taxonmap reference_proteomes.taxid_map --taxonnodes /scicore/home/ebertd/angpas00/hifi/new_taxdump/nodes.dmp --taxonnames /scicore/home/ebertd/angpas00/hifi/new_taxdump/names.dmp -d reference_proteomes.dmnd

diamond blastx \
      --query longi.p_ctg.fa \
      --db ../uniprot/reference_proteomes.dmnd \
      --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
      --sensitive \
      --max-target-seqs 1 \
      --evalue 1e-25 \
      --threads 50 \
      -F 15 -b4 -c1 \
      > diamond-scicore.out  

rm -r ../uniprot/
        
# map reads to assembly to get coverage per contig

minimap2 -ax map-hifi -t 54 longi.p_ctg.fa ../Pool1_Circular_Consensus_Sequencing_Reads/m64156_230621_203451.bc2096--bc2096.hifi_reads.fastq.gz > full.sam

samtools view -Sb -@ 10 full.sam | samtools sort -@ 10 -T abcxyz -o full.bam
rm full.sam
samtools index -b full.bam       
                
# 4) summarize info with blobtools        
        
# run blobtools
blobtools nodesdb --nodes /home/pascal/assemblies/denti/new_taxdump/nodes.dmp --names /home/pascal/assemblies/denti/new_taxdump/names.dmp

blobtools create -i longi.p_ctg.fa -b full.bam \
 -t diamond-scicore.out -t blast.scicore.out -o FULL.blobtools \
 --nodes /home/pascal/assemblies/denti/new_taxdump/nodes.dmp --names /home/pascal/assemblies/denti/new_taxdump/names.dmp

blobtools view -i FULL.blobtools.blobDB.json
blobtools plot -i FULL.blobtools.blobDB.json




# filter based on visual and numerical inspection:

grep -v [0-9]c FULL.blobtools.blobDB.bestsum.table.txt | grep -v "#" | awk '{if ($3 > 0.35 && $3 < 0.41 && $5 > 11) print $0}' | grep -v Pseudomonadota | cut -f1 > denti.names 
grep -v [0-9]c FULL.blobtools.blobDB.bestsum.table.txt | grep -v "#" | awk '{if ($3 > 0.35 && $3 < 0.41 && $5 > 15) print $0}' | grep -v Pseudomonadota | cut -f1 > longi.names
 



# repeat blobtools analysis with updated assembly

grep -f longi.names blast.scicore.out > blast.longi.out
grep -f longi.names diamond-scicore.out > diamond.longi.out

minimap2 -ax map-hifi -t 54 longi.fa ../Pool1_Circular_Consensus_Sequencing_Reads/m64156_230621_203451.bc2096--bc2096.hifi_reads.fastq.gz > longi.sam

samtools view -Sb -@ 10 longi.sam | samtools sort -@ 10 -T abcxyz -o longi.bam
rm longi.sam
samtools index -b longi.bam

blobtools create -i longi.fa -b longi.bam \
 -t diamond.longi.out -t blast.longi.out -o FULL.blobtools \
 --nodes /home/pascal/assemblies/denti/new_taxdump/nodes.dmp --names /home/pascal/assemblies/denti/new_taxdump/names.dmp

blobtools view -i FULL.blobtools.blobDB.json
blobtools plot -i FULL.blobtools.blobDB.json



# look for mtDNA

denti: ptg000049c	16550	0.3595	0	286.6361	Arthropoda	790.0	0
longi: no circular molecule identified as mtDNA