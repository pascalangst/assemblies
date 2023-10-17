# D. magna from island Spicarna (SP1) infected with Spirobacillus

# collected on the 11.09.2023
# 4 bottles = 4 populations pooled for ONT sequencing (bottles Nr. 21, 40, 202, 35)
# ~ 30 animals in total
# puregene extraction protocol without glycogen but with mesh for pipetting supernatant after protK digest

# base calling of ONT reads (dorado in duplex mode)
singularity run --nv  /export/soft/singularity-containers/dorado/dorado_0.3.4.sif dorado duplex /scicore/home/ebertd/angpas00/software/dorado/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 -v pod5_pass/ > duplex_231012.bam

# bam to fq
module load BEDTools/2.30.0-GCC-10.3.0 
bedtools bamtofastq -i duplex_231012.bam -fq duplex_231012.fq

# assembly (flye in meta mode)
flye --nano-raw duplex_231012.fq.gz --meta --threads 55 --out-dir flye_meta_ont

# extract single contig which represents circular chromosome of Spirobacillus and annotate it

# set up database for bakta
#wget https://zenodo.org/record/7669534/files/db.tar.gz
#tar -xzf db.tar.gz
#rm db.tar.gz
#amrfinder_update --force_update --database db/amrfinderplus-db/

# annotation using bakta
bakta --db ~/bioinformatics/bakta/db --output bakta/ --prefix Spiro2 --locus-tag Spiro2 --threads 8 assembly_meta_ont.Spiro.rotated.fasta 

#parse genome sequences...
#	imported: 1
#	filtered & revised: 1
#	contigs: 1
#
#start annotation...
#predict tRNAs...
#	found: 42
#predict tmRNAs...
#	found: 1
#predict rRNAs...
#	found: 13
#predict ncRNAs...
#	found: 6
#predict ncRNA regions...
#	found: 0
#predict CRISPR arrays...
#	found: 2
#predict & annotate CDSs...
#	predicted: 2499 
#	discarded spurious: 5
#	revised translational exceptions: 0
#	detected IPSs: 683
#	found PSCs: 1308
#	found PSCCs: 184
#	lookup annotations...
#	conduct expert systems...
#		amrfinder: 0
#		protein sequences: 0
#	combine annotations and mark hypotheticals...
#	detect pseudogenes...
#		pseudogene candidates: 48
#		found pseudogenes: 38
#analyze hypothetical proteins: 320
#	detected Pfam hits: 66 
#	calculated proteins statistics
#	revise special cases...
#extract sORF...
#	potential: 40093
#	discarded due to overlaps: 35757
#	discarded spurious: 0
#	detected IPSs: 0
#	found PSCs: 0
#	lookup annotations...
#	filter and combine annotations...
#	filtered sORFs: 0
#detect gaps...
#	found: 0
#detect oriCs/oriVs...
#	found: 0
#detect oriTs...
#	found: 0
#apply feature overlap filters...
#select features and create locus tags...
#selected: 2550
#improve annotations...
#	revised gene symbols: 0
#
#genome statistics:
#	Genome size: 2,806,830 bp
#	Contigs/replicons: 1
#	GC: 32.2 %
#	N50: 2,806,830
#	N ratio: 0.0 %
#	coding density: 91.1 %
#
#annotation summary:
#	tRNAs: 42
#	tmRNAs: 1
#	rRNAs: 13
#	ncRNAs: 6
#	ncRNA regions: 0
#	CRISPR arrays: 2
#	CDSs: 2486
#		hypotheticals: 320
#		pseudogenes: 38
#		signal peptides: 0
#	sORFs: 0
#	gaps: 0
#	oriCs/oriVs: 0
#	oriTs: 0
#
#export annotation results to: /home/pascal/assemblies/Spiro/ont/bakta
#	human readable TSV...
#	GFF3...
#	INSDC GenBank & EMBL...
#	genome sequences...
#	feature nucleotide sequences...
#	translated CDS sequences...
#	circular genome plot...
#	hypothetical TSV...
#	translated hypothetical CDS sequences...
#	machine readable JSON...
#	genome and annotation summary...
#
#If you use these results please cite Bakta: https://doi.org/10.1099/mgen.0.000685
#Annotation successfully finished in 9:33 [mm:ss].