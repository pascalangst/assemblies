### Comparative Genomics with newly generated Assemblies

# first, assemblies need to be annotated

## Binucleata daphniae

cd /home/pascal/BE-OM-3_Bd-2-2_depl/CA.mr.49.17.15.0.02/

conda activate mamba-funannotate-Ht # 1.8.14

funannotate clean -i curated2.fasta -o BE-OM-3.v1.clean.fa
funannotate sort -i BE-OM-3.v1.clean.fa -o BE-OM-3.v1.sort.fa
funannotate mask -i BE-OM-3.v1.sort.fa -o BE-OM-3.v1.mask.fa



export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate-Ht/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate-Ht
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate-Ht/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate-Ht/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate-Ht/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate-Ht/opt/pasa-2.5.2

# run funannotate with fungal gene hints
funannotate predict -i BE-OM-3.v1.mask.fa -o funannotate_BE-OM-3.v1 -s "Binucleata daphniae" --isolate BE-OM-3 --augustus_species encephalitozoon_cuniculi_GB --busco_db microsporidia --busco_seed_species encephalitozoon_cuniculi_GB --cpus 12 \
--transcript_evidence \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/GCA_000442015.1_Rozella_k41_t100_rna_from_genomic.fna.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_AalgeraePRA109_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_AalgeraePRA339_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_EaedisUSNM41457_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_EcuniculiGBM1_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_MdaphniaeUGP3_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NausubeliERTm2_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NausubeliERTm6_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NbombycisCQ1_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NceranaeBRL01_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NparisiiERTm1_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NparisiiERTm3_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NparisiiERTm3_AnnotatedTranscripts.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_OcolligataOC4_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_PneurophiliaMK1_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_Slophii42_110_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_ThominisUnknown_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_VcorneaeATCC50505_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_Vculicisfloridensis_AnnotatedTranscripts.fasta.oneline.fasta

#-------------------------------------------------------
#[Mar 23 08:34 AM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.13
#[Mar 23 08:34 AM]: Running funannotate v1.8.14
#[Mar 23 08:34 AM]: Skipping CodingQuarry as no --rna_bam passed
#[Mar 23 08:34 AM]: Parsed training data, run ab-initio gene predictors as follows:
#  Program      Training-Method
#  augustus     pretrained     
#  genemark     selftraining   
#  glimmerhmm   busco          
#  snap         busco          
#[Mar 23 08:34 AM]: Loading genome assembly and parsing soft-masked repetitive sequences
#[Mar 23 08:34 AM]: Genome loaded: 31 scaffolds; 9,578,088 bp; 9.70% repeats masked
#[Mar 23 08:34 AM]: Aligning transcript evidence to genome with minimap2
#[Mar 23 08:34 AM]: Found 5 alignments, wrote GFF3 and Augustus hints to file
#[Mar 23 08:34 AM]: Mapping 555,555 proteins to genome using diamond and exonerate
#[Mar 23 08:36 AM]: Found 108,475 preliminary alignments with diamond in 0:00:49 --> generated FASTA files for exonerate in 0:00:17
#     Progress: 108475 complete, 0 failed, 0 remaining           
#[Mar 23 08:45 AM]: Exonerate finished in 0:09:49: found 314 alignments
#[Mar 23 08:46 AM]: Running GeneMark-ES on assembly
#[Mar 23 08:50 AM]: 4,291 predictions from GeneMark
#[Mar 23 08:50 AM]: Running BUSCO to find conserved gene models for training ab-initio predictors
#[Mar 23 08:53 AM]: 184 valid BUSCO predictions found, validating protein sequences
#[Mar 23 08:54 AM]: 184 BUSCO predictions validated
#[Mar 23 08:54 AM]: Running Augustus gene prediction using encephalitozoon_cuniculi_GB parameters
#     Progress: 35 complete, 0 failed, 0 remaining         
#[Mar 23 08:55 AM]: 35 predictions from Augustus
#[Mar 23 08:55 AM]: Pulling out high quality Augustus predictions
#[Mar 23 08:55 AM]: Found 3 high quality predictions from Augustus (>90% exon evidence)
#[Mar 23 08:55 AM]: Running SNAP gene prediction, using training data: funannotate_BE-OM-3.v1/predict_misc/busco.final.gff3
#[Mar 23 08:56 AM]: 3,896 predictions from SNAP
#[Mar 23 08:56 AM]: Running GlimmerHMM gene prediction, using training data: funannotate_BE-OM-3.v1/predict_misc/busco.final.gff3
#[Mar 23 08:57 AM]: 3,641 predictions from GlimmerHMM
#[Mar 23 08:57 AM]: Summary of gene models passed to EVM (weights):
#  Source         Weight   Count
#  Augustus       1        32   
#  Augustus HiQ   2        3    
#  GeneMark       1        4291 
#  GlimmerHMM     1        3641 
#  snap           1        3896 
#  Total          -        11863
#[Mar 23 08:57 AM]: EVM: partitioning input to ~ 35 genes per partition using min 1500 bp interval
#     Progress: 120 complete, 0 failed, 0 remaining         
#[Mar 23 09:12 AM]: Converting to GFF3 and collecting all EVM results
#[Mar 23 09:12 AM]: 3,333 total gene models from EVM
#[Mar 23 09:12 AM]: Generating protein fasta files from 3,333 EVM models
#[Mar 23 09:13 AM]: now filtering out bad gene models (< 50 aa in length, transposable elements, etc).
#[Mar 23 09:13 AM]: Found 486 gene models to remove: 0 too short; 0 span gaps; 486 transposable elements
#[Mar 23 09:13 AM]: 2,847 gene models remaining
#[Mar 23 09:13 AM]: Predicting tRNAs
#[Mar 23 09:13 AM]: 61 tRNAscan models are valid (non-overlapping)
#[Mar 23 09:13 AM]: Generating GenBank tbl annotation file
#[Mar 23 09:13 AM]: Collecting final annotation files for 2,908 total gene models
#[Mar 23 09:13 AM]: Converting to final Genbank format
#[Mar 23 09:14 AM]: Funannotate predict is finished, output files are in the funannotate_BE-OM-3.v1/predict_results folder
#[Mar 23 09:14 AM]: Your next step might be functional annotation, suggested commands:
#-------------------------------------------------------
#Run InterProScan (manual install): 
#funannotate iprscan -i funannotate_BE-OM-3.v1 -c 12
#
#Run antiSMASH (optional): 
#funannotate remote -i funannotate_BE-OM-3.v1 -m antismash -e youremail@server.edu
#
#Annotate Genome: 
#funannotate annotate -i funannotate_BE-OM-3.v1 --cpus 12 --sbt yourSBTfile.txt
#-------------------------------------------------------
#                
#[Mar 23 09:14 AM]: Training parameters file saved: funannotate_BE-OM-3.v1/predict_results/encephalitozoon_cuniculi_GB.parameters.json
#[Mar 23 09:14 AM]: Add species parameters to database:
#
#  funannotate species -s encephalitozoon_cuniculi_GB -a funannotate_BE-OM-3.v1/predict_results/encephalitozoon_cuniculi_GB.parameters.json
#


mkdir funannotate_BE-OM-3.v1/annotate_misc/
conda activate interproscan
python ~/Ht_methylation_project/analysis/iso-seq/long_read_protocol/iprscan-local.py -i funannotate_BE-OM-3.v1/ -m local -c 12

conda deactivate
export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/pasa-2.4.1

funannotate annotate -i funannotate_BE-OM-3.v1/ --cpus 12 --iprscan funannotate_BE-OM-3.v1/annotate_misc/iprscan.xml --busco_db microsporidia

#-------------------------------------------------------
#[Mar 23 09:49 AM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.13
#[Mar 23 09:49 AM]: Running 1.8.14
#[Mar 23 09:49 AM]: No NCBI SBT file given, will use default, however if you plan to submit to NCBI, create one and pass it here '--sbt'
#[Mar 23 09:49 AM]: Parsing input files
#[Mar 23 09:49 AM]: Existing tbl found: funannotate_BE-OM-3.v1/predict_results/Binucleata_daphniae_BE-OM-3.tbl
#[Mar 23 09:49 AM]: Adding Functional Annotation to Binucleata daphniae, NCBI accession: None
#[Mar 23 09:49 AM]: Annotation consists of: 2,908 gene models
#[Mar 23 09:49 AM]: 2,847 protein records loaded
#[Mar 23 09:49 AM]: Running HMMer search of PFAM version 35.0
#[Mar 23 09:51 AM]: 2,065 annotations added
#[Mar 23 09:51 AM]: Running Diamond blastp search of UniProt DB version 2022_01
#[Mar 23 09:51 AM]: 45 valid gene/product annotations from 55 total
#[Mar 23 09:51 AM]: Running Eggnog-mapper
#[Mar 23 10:03 AM]: Parsing EggNog Annotations
#[Mar 23 10:03 AM]: EggNog version parsed as 2.1.9
#[Mar 23 10:03 AM]: 3,298 COG and EggNog annotations added
#[Mar 23 10:03 AM]: Combining UniProt/EggNog gene and product names using Gene2Product version 1.76
#[Mar 23 10:03 AM]: 805 gene name and product description annotations added
#[Mar 23 10:03 AM]: Running Diamond blastp search of MEROPS version 12.0
#[Mar 23 10:03 AM]: 64 annotations added
#[Mar 23 10:03 AM]: Annotating CAZYmes using HMMer search of dbCAN version 10.0
#[Mar 23 10:03 AM]: 11 annotations added
#[Mar 23 10:03 AM]: Annotating proteins with BUSCO microsporidia models
#[Mar 23 10:03 AM]: 484 annotations added
#[Mar 23 10:03 AM]: Predicting secreted and transmembrane proteins using Phobius
#     Progress: 2847 complete, 0 failed, 0 remaining           
#[Mar 23 10:04 AM]: Predicting secreted proteins with SignalP
#[Mar 23 10:10 AM]: 111 secretome and 357 transmembane annotations added
#[Mar 23 10:10 AM]: Parsing InterProScan5 XML file
#[Mar 23 10:10 AM]: Found 0 duplicated annotations, adding 8,149 valid annotations
#[Mar 23 10:10 AM]: Converting to final Genbank format, good luck!
#[Mar 23 10:10 AM]: Creating AGP file and corresponding contigs file
#[Mar 23 10:10 AM]: Writing genome annotation table.
#[Mar 23 10:10 AM]: Funannotate annotate has completed successfully!
#
#        We need YOUR help to improve gene names/product descriptions:
#           0 gene/products names MUST be fixed, see funannotate_BE-OM-3.v1/annotate_results/Gene2Products.must-fix.txt
#           0 gene/product names need to be curated, see funannotate_BE-OM-3.v1/annotate_results/Gene2Products.need-curating.txt
#           18 gene/product names passed but are not in Database, see funannotate_BE-OM-3.v1/annotate_results/Gene2Products.new-names-passed.txt
#
#        Please consider contributing a PR at https://github.com/nextgenusfs/gene2product
#
#-------------------------------------------------------














## Gurleya daphniae

cd /home/pascal/Gv_SP/nextDenovo_polished_sgs-lgs_medaka/



funannotate clean -i curated_without_ctg000290_np1212_pilon.fasta  -o Gv_SP.v2.clean.fa
funannotate sort -i Gv_SP.v2.clean.fa -o Gv_SP.v2.sort.fa
funannotate mask -i Gv_SP.v2.sort.fa -o Gv_SP.v2.mask.fa



export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate-Ht/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate-Ht
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate-Ht/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate-Ht/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate-Ht/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate-Ht/opt/pasa-2.5.2

# run funannotate with fungal gene hints
funannotate predict -i Gv_SP.v2.mask.fa -o funannotate_Gv_SP.v2 -s "Gurleya daphniae" --isolate Gv_SP --augustus_species encephalitozoon_cuniculi_GB --busco_db microsporidia --busco_seed_species encephalitozoon_cuniculi_GB --cpus 12 \
--transcript_evidence \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/GCA_000442015.1_Rozella_k41_t100_rna_from_genomic.fna.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_AalgeraePRA109_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_AalgeraePRA339_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_EaedisUSNM41457_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_EcuniculiGBM1_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_MdaphniaeUGP3_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NausubeliERTm2_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NausubeliERTm6_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NbombycisCQ1_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NceranaeBRL01_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NparisiiERTm1_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NparisiiERTm3_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NparisiiERTm3_AnnotatedTranscripts.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_OcolligataOC4_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_PneurophiliaMK1_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_Slophii42_110_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_ThominisUnknown_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_VcorneaeATCC50505_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_Vculicisfloridensis_AnnotatedTranscripts.fasta.oneline.fasta

#-------------------------------------------------------
#[Jun 13 08:46 AM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.13
#[Jun 13 08:46 AM]: Running funannotate v1.8.14
#
#[Jun 13 08:46 AM]: Skipping CodingQuarry as no --rna_bam passed
#[Jun 13 08:46 AM]: Parsed training data, run ab-initio gene predictors as follows:
#  Program      Training-Method
#  augustus     pretrained
#  genemark     selftraining
#  glimmerhmm   busco
#  snap         busco
#[Jun 13 08:46 AM]: Loading genome assembly and parsing soft-masked repetitive sequences
#[Jun 13 08:46 AM]: Genome loaded: 20 scaffolds; 16,375,077 bp; 22.21% repeats masked
#[Jun 13 08:47 AM]: Aligning transcript evidence to genome with minimap2
#[Jun 13 08:47 AM]: Found 9 alignments, wrote GFF3 and Augustus hints to file
#[Jun 13 08:47 AM]: Mapping 555,555 proteins to genome using diamond and exonerate
#[Jun 13 08:48 AM]: Found 109,130 preliminary alignments with diamond in 0:00:54 --> generated FASTA files for exonerate in 0:00:17
#     Progress: 109130 complete, 0 failed, 0 remaining
#[Jun 13 08:58 AM]: Exonerate finished in 0:09:54: found 358 alignments
#[Jun 13 08:58 AM]: Running GeneMark-ES on assembly
#[Jun 13 09:13 AM]: 5,379 predictions from GeneMark
#[Jun 13 09:13 AM]: Running BUSCO to find conserved gene models for training ab-initio predictors
#[Jun 13 09:16 AM]: 223 valid BUSCO predictions found, validating protein sequences
#[Jun 13 09:16 AM]: 223 BUSCO predictions validated
#[Jun 13 09:16 AM]: Running Augustus gene prediction using encephalitozoon_cuniculi_GB parameters
#     Progress: 42 complete, 0 failed, 0 remaining
#[Jun 13 09:18 AM]: 10 predictions from Augustus
#[Jun 13 09:18 AM]: Pulling out high quality Augustus predictions
#[Jun 13 09:18 AM]: Found 3 high quality predictions from Augustus (>90% exon evidence)
#[Jun 13 09:18 AM]: Running SNAP gene prediction, using training data: funannotate_Gv_SP.v2/predict_misc/busco.final.gff3
#[Jun 13 09:19 AM]: 4,328 predictions from SNAP
#[Jun 13 09:19 AM]: Running GlimmerHMM gene prediction, using training data: funannotate_Gv_SP.v2/predict_misc/busco.final.gff3
#[Jun 13 09:21 AM]: 4,709 predictions from GlimmerHMM
#[Jun 13 09:21 AM]: Summary of gene models passed to EVM (weights):
#  Source         Weight   Count
#  Augustus       1        7
#  Augustus HiQ   2        3
#  GeneMark       1        5379
#  GlimmerHMM     1        4709
#  snap           1        4328
#  Total          -        14426
#[Jun 13 09:21 AM]: EVM: partitioning input to ~ 35 genes per partition using min 1500 bp interval
#     Progress: 165 complete, 0 failed, 0 remaining
#[Jun 13 09:28 AM]: Converting to GFF3 and collecting all EVM results
#[Jun 13 09:28 AM]: 3,976 total gene models from EVM
#[Jun 13 09:28 AM]: Generating protein fasta files from 3,976 EVM models
#[Jun 13 09:28 AM]: now filtering out bad gene models (< 50 aa in length, transposable elements, etc).
#[Jun 13 09:28 AM]: Found 433 gene models to remove: 0 too short; 0 span gaps; 433 transposable elements
#[Jun 13 09:28 AM]: 3,543 gene models remaining
#[Jun 13 09:28 AM]: Predicting tRNAs
#[Jun 13 09:28 AM]: 62 tRNAscan models are valid (non-overlapping)
#[Jun 13 09:28 AM]: Generating GenBank tbl annotation file
#[Jun 13 09:28 AM]: Collecting final annotation files for 3,605 total gene models
#[Jun 13 09:28 AM]: Converting to final Genbank format
#[Jun 13 09:29 AM]: Funannotate predict is finished, output files are in the funannotate_Gv_SP.v2/predict_results folder
#[Jun 13 09:29 AM]: Your next step might be functional annotation, suggested commands:
#-------------------------------------------------------
#Run InterProScan (manual install):
#funannotate iprscan -i funannotate_Gv_SP.v2 -c 12
#
#Run antiSMASH (optional):
#funannotate remote -i funannotate_Gv_SP.v2 -m antismash -e youremail@server.edu
#
#Annotate Genome:
#funannotate annotate -i funannotate_Gv_SP.v2 --cpus 12 --sbt yourSBTfile.txt
#-------------------------------------------------------
#
#[Jun 13 09:29 AM]: Training parameters file saved: funannotate_Gv_SP.v2/predict_results/encephalitozoon_cuniculi_GB.parameters.json
#[Jun 13 09:29 AM]: Add species parameters to database:
#
#  funannotate species -s encephalitozoon_cuniculi_GB -a funannotate_Gv_SP.v2/predict_results/encephalitozoon_cuniculi_GB.parameters.json

mkdir funannotate_Gv_SP.v2/annotate_misc/
conda activate interproscan
python ~/Ht_methylation_project/analysis/iso-seq/long_read_protocol/iprscan-local.py -i funannotate_Gv_SP.v2/ -m local -c 12

conda deactivate
export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/pasa-2.4.1

funannotate annotate -i funannotate_Gv_SP.v2/ --cpus 12 --iprscan funannotate_Gv_SP.v2/annotate_misc/iprscan.xml --busco_db microsporidia
#-------------------------------------------------------
#[Jun 13 10:08 AM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.13
#[Jun 13 10:08 AM]: Running 1.8.14
#[Jun 13 10:08 AM]: No NCBI SBT file given, will use default, however if you plan to submit to NCBI, create one and pass it here '--sbt'
#[Jun 13 10:08 AM]: Parsing input files
#[Jun 13 10:08 AM]: Existing tbl found: funannotate_Gv_SP.v2/predict_results/Gurleya_daphniae_Gv_SP.tbl
#[Jun 13 10:08 AM]: Adding Functional Annotation to Gurleya daphniae, NCBI accession: None
#[Jun 13 10:08 AM]: Annotation consists of: 3,605 gene models
#[Jun 13 10:08 AM]: 3,543 protein records loaded
#[Jun 13 10:08 AM]: Running HMMer search of PFAM version 35.0
#[Jun 13 10:09 AM]: 1,945 annotations added
#[Jun 13 10:09 AM]: Running Diamond blastp search of UniProt DB version 2022_01
#[Jun 13 10:10 AM]: 44 valid gene/product annotations from 60 total
#[Jun 13 10:10 AM]: Running Eggnog-mapper
#[Jun 13 10:22 AM]: Parsing EggNog Annotations
#[Jun 13 10:22 AM]: EggNog version parsed as 2.1.9
#[Jun 13 10:22 AM]: 3,347 COG and EggNog annotations added
#[Jun 13 10:22 AM]: Combining UniProt/EggNog gene and product names using Gene2Product version 1.76
#[Jun 13 10:22 AM]: 782 gene name and product description annotations added
#[Jun 13 10:22 AM]: Running Diamond blastp search of MEROPS version 12.0
#[Jun 13 10:22 AM]: 71 annotations added
#[Jun 13 10:22 AM]: Annotating CAZYmes using HMMer search of dbCAN version 10.0
#[Jun 13 10:23 AM]: 9 annotations added
#[Jun 13 10:23 AM]: Annotating proteins with BUSCO microsporidia models
#[Jun 13 10:23 AM]: 463 annotations added
#[Jun 13 10:23 AM]: Predicting secreted and transmembrane proteins using Phobius
#     Progress: 3543 complete, 0 failed, 0 remaining
#[Jun 13 10:23 AM]: Predicting secreted proteins with SignalP
#[Jun 13 10:31 AM]: 158 secretome and 535 transmembane annotations added
#[Jun 13 10:31 AM]: Parsing InterProScan5 XML file
#[Jun 13 10:31 AM]: Found 0 duplicated annotations, adding 8,236 valid annotations
#[Jun 13 10:31 AM]: Converting to final Genbank format, good luck!
#[Jun 13 10:32 AM]: Creating AGP file and corresponding contigs file
#[Jun 13 10:32 AM]: Writing genome annotation table.
#[Jun 13 10:32 AM]: Funannotate annotate has completed successfully!
#
#        We need YOUR help to improve gene names/product descriptions:
#           0 gene/products names MUST be fixed, see funannotate_Gv_SP.v2/annotate_results/Gene2Products.must-fix.txt
#           2 gene/product names need to be curated, see funannotate_Gv_SP.v2/annotate_results/Gene2Products.need-curating.txt
#           17 gene/product names passed but are not in Database, see funannotate_Gv_SP.v2/annotate_results/Gene2Products.new-names-passed.txt
#
#        Please consider contributing a PR at https://github.com/nextgenusfs/gene2product
#
#-------------------------------------------------------










## Larssonia obtusa
# changed genemark parameters for selftraining (use contigs > 10kb instead > 50kb)
# dirty solution: cnage python code in lib/python3.8/site-packages/funannotate/library.py, function RunGeneMarkES, to "command, '--ES', '--min_contig', str(10000), '--max_intron', str(maxintron), ..."
# don't forget to remove '--min_contig', str(10000), after use!

cd /home/pascal/assemblies/SK-39_Larssonia/megahit/



funannotate clean -i contigs.cov50.500.GC33.noSMR.fasta -o FI-SK-39.v1.clean.fa
funannotate sort -i FI-SK-39.v1.clean.fa -o FI-SK-39.v1.sort.fa
funannotate mask -i FI-SK-39.v1.sort.fa -o FI-SK-39.v1.mask.fa



export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate-Ht/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate-Ht
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate-Ht/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate-Ht/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate-Ht/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate-Ht/opt/pasa-2.5.2

# run funannotate with fungal gene hints
funannotate predict -i FI-SK-39.v1.mask.fa -o funannotate_FI-SK-39.v1.1 -s "Larssonia obtusa" --isolate FI-SK-39 --augustus_species encephalitozoon_cuniculi_GB --busco_db microsporidia --busco_seed_species encephalitozoon_cuniculi_GB --cpus 12 \
--transcript_evidence \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/GCA_000442015.1_Rozella_k41_t100_rna_from_genomic.fna.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_AalgeraePRA109_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_AalgeraePRA339_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_EaedisUSNM41457_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_EcuniculiGBM1_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_MdaphniaeUGP3_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NausubeliERTm2_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NausubeliERTm6_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NbombycisCQ1_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NceranaeBRL01_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NparisiiERTm1_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NparisiiERTm3_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_NparisiiERTm3_AnnotatedTranscripts.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_OcolligataOC4_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_PneurophiliaMK1_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_Slophii42_110_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_ThominisUnknown_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_VcorneaeATCC50505_AnnotatedTranscripts.fasta.oneline.fasta \
/home/pascal/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/MicrosporidiaDB-56_Vculicisfloridensis_AnnotatedTranscripts.fasta.oneline.fasta
#-------------------------------------------------------
#[Jun 13 11:15 AM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.13
#[Jun 13 11:15 AM]: Running funannotate v1.8.14
#[Jun 13 11:15 AM]: Skipping CodingQuarry as no --rna_bam passed
#[Jun 13 11:15 AM]: Parsed training data, run ab-initio gene predictors as follows:
#  Program      Training-Method
#  augustus     pretrained
#  genemark     selftraining
#  glimmerhmm   busco
#  snap         busco
#[Jun 13 11:15 AM]: Loading genome assembly and parsing soft-masked repetitive sequences
#[Jun 13 11:15 AM]: Genome loaded: 5,654 scaffolds; 18,130,460 bp; 11.08% repeats masked
#[Jun 13 11:15 AM]: Existing transcript alignments found: funannotate_FI-SK-39.v1.1/predict_misc/transcript_alignments.gff3
#[Jun 13 11:15 AM]: Existing protein alignments found: funannotate_FI-SK-39.v1.1/predict_misc/protein_alignments.gff3
#[Jun 13 11:15 AM]: Running GeneMark-ES on assembly
#[Jun 13 11:34 AM]: 6,328 predictions from GeneMark
#[Jun 13 11:34 AM]: Running BUSCO to find conserved gene models for training ab-initio predictors
#[Jun 13 11:36 AM]: 227 valid BUSCO predictions found, validating protein sequences
#[Jun 13 11:37 AM]: 226 BUSCO predictions validated
#[Jun 13 11:37 AM]: Running Augustus gene prediction using encephalitozoon_cuniculi_GB parameters
#     Progress: 5654 complete, 0 failed, 0 remaining
#[Jun 13 11:43 AM]: 33 predictions from Augustus
#[Jun 13 11:43 AM]: Pulling out high quality Augustus predictions
#[Jun 13 11:43 AM]: Found 1 high quality predictions from Augustus (>90% exon evidence)
#[Jun 13 11:43 AM]: Running SNAP gene prediction, using training data: funannotate_FI-SK-39.v1.1/predict_misc/busco.final.gff3
#[Jun 13 11:44 AM]: 6,426 predictions from SNAP
#[Jun 13 11:44 AM]: Running GlimmerHMM gene prediction, using training data: funannotate_FI-SK-39.v1.1/predict_misc/busco.final.gff3
#[Jun 13 11:47 AM]: 4,644 predictions from GlimmerHMM
#[Jun 13 11:47 AM]: Summary of gene models passed to EVM (weights):
#  Source         Weight   Count
#  Augustus       1        32
#  Augustus HiQ   2        1
#  GeneMark       1        6328
#  GlimmerHMM     1        4644
#  snap           1        6426
#  Total          -        17431
#[Jun 13 11:47 AM]: EVM: partitioning input to ~ 35 genes per partition using min 1500 bp interval
#     Progress: 3772 complete, 0 failed, 0 remaining
#[Jun 13 11:54 AM]: Converting to GFF3 and collecting all EVM results
#[Jun 13 11:54 AM]: 4,876 total gene models from EVM
#[Jun 13 11:54 AM]: Generating protein fasta files from 4,876 EVM models
#[Jun 13 11:54 AM]: now filtering out bad gene models (< 50 aa in length, transposable elements, etc).
#[Jun 13 11:54 AM]: Found 645 gene models to remove: 8 too short; 0 span gaps; 637 transposable elements
#[Jun 13 11:54 AM]: 4,231 gene models remaining
#[Jun 13 11:54 AM]: Predicting tRNAs
#[Jun 13 12:07 PM]: 60 tRNAscan models are valid (non-overlapping)
#[Jun 13 12:07 PM]: Generating GenBank tbl annotation file
#[Jun 13 12:07 PM]: Collecting final annotation files for 4,291 total gene models
#[Jun 13 12:07 PM]: Converting to final Genbank format
#[Jun 13 12:08 PM]: Funannotate predict is finished, output files are in the funannotate_FI-SK-39.v1.1/predict_results folder
#[Jun 13 12:08 PM]: Your next step might be functional annotation, suggested commands:
#-------------------------------------------------------
#Run InterProScan (manual install):
#funannotate iprscan -i funannotate_FI-SK-39.v1.1 -c 12
#
#Run antiSMASH (optional):
#funannotate remote -i funannotate_FI-SK-39.v1.1 -m antismash -e youremail@server.edu
#
#Annotate Genome:
#funannotate annotate -i funannotate_FI-SK-39.v1.1 --cpus 12 --sbt yourSBTfile.txt
#-------------------------------------------------------
#
#[Jun 13 12:08 PM]: Training parameters file saved: funannotate_FI-SK-39.v1.1/predict_results/encephalitozoon_cuniculi_GB.parameters.json
#[Jun 13 12:08 PM]: Add species parameters to database:
#
#  funannotate species -s encephalitozoon_cuniculi_GB -a funannotate_FI-SK-39.v1.1/predict_results/encephalitozoon_cuniculi_GB.parameters.json

mkdir funannotate_FI-SK-39.v1.1/annotate_misc/
conda activate interproscan
python ~/Ht_methylation_project/analysis/iso-seq/long_read_protocol/iprscan-local.py -i funannotate_FI-SK-39.v1.1/ -m local -c 12

conda deactivate
export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/pasa-2.4.1

funannotate annotate -i funannotate_FI-SK-39.v1.1/ --cpus 12 --iprscan funannotate_FI-SK-39.v1.1/annotate_misc/iprscan.xml --busco_db microsporidia
#-------------------------------------------------------
#[Jun 13 12:20 PM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.13
#[Jun 13 12:20 PM]: Running 1.8.14
#[Jun 13 12:20 PM]: No NCBI SBT file given, will use default, however if you plan to submit to NCBI, create one and pass it here '--sbt'
#[Jun 13 12:20 PM]: Parsing input files
#[Jun 13 12:20 PM]: Existing tbl found: funannotate_FI-SK-39.v1.1/predict_results/Larssonia_obtusa_FI-SK-39.tbl
#[Jun 13 12:20 PM]: Adding Functional Annotation to Conglomerata obtusa, NCBI accession: None
#[Jun 13 12:20 PM]: Annotation consists of: 4,291 gene models
#[Jun 13 12:20 PM]: 4,231 protein records loaded
#[Jun 13 12:20 PM]: Running HMMer search of PFAM version 35.0
#[Jun 13 12:22 PM]: 1,914 annotations added
#[Jun 13 12:22 PM]: Running Diamond blastp search of UniProt DB version 2022_01
#[Jun 13 12:23 PM]: 40 valid gene/product annotations from 52 total
#[Jun 13 12:23 PM]: Running Eggnog-mapper
#[Jun 13 12:34 PM]: Parsing EggNog Annotations
#[Jun 13 12:34 PM]: EggNog version parsed as 2.1.9
#[Jun 13 12:34 PM]: 3,432 COG and EggNog annotations added
#[Jun 13 12:34 PM]: Combining UniProt/EggNog gene and product names using Gene2Product version 1.76
#[Jun 13 12:34 PM]: 746 gene name and product description annotations added
#[Jun 13 12:34 PM]: Running Diamond blastp search of MEROPS version 12.0
#[Jun 13 12:34 PM]: 67 annotations added
#[Jun 13 12:34 PM]: Annotating CAZYmes using HMMer search of dbCAN version 10.0
#[Jun 13 12:34 PM]: 10 annotations added
#[Jun 13 12:34 PM]: Annotating proteins with BUSCO microsporidia models
#[Jun 13 12:34 PM]: 418 annotations added
#[Jun 13 12:34 PM]: Predicting secreted and transmembrane proteins using Phobius
#     Progress: 4231 complete, 0 failed, 0 remaining
#[Jun 13 12:35 PM]: Predicting secreted proteins with SignalP
#[Jun 13 12:44 PM]: 175 secretome and 556 transmembane annotations added
#[Jun 13 01:35 PM]: Parsing InterProScan5 XML file
#[Jun 13 01:35 PM]: Found 0 duplicated annotations, adding 8,177 valid annotations
#[Jun 13 01:35 PM]: Converting to final Genbank format, good luck!
#[Jun 13 01:36 PM]: Creating AGP file and corresponding contigs file
#[Jun 13 01:36 PM]: Writing genome annotation table.
#[Jun 13 01:36 PM]: Funannotate annotate has completed successfully!
#
#        We need YOUR help to improve gene names/product descriptions:
#           0 gene/products names MUST be fixed, see funannotate_FI-SK-39.v1.1/annotate_results/Gene2Products.must-fix.txt
#           2 gene/product names need to be curated, see funannotate_FI-SK-39.v1.1/annotate_results/Gene2Products.need-curating.txt
#           22 gene/product names passed but are not in Database, see funannotate_FI-SK-39.v1.1/annotate_results/Gene2Products.new-names-passed.txt
#
#        Please consider contributing a PR at https://github.com/nextgenusfs/gene2product
#
#-------------------------------------------------------


















## comparative analysis

for file in ~/assemblies/Mito_SP/masurca_ont_pb/CA.mr.99.17.15.0.02/contigs.SP.braker.dikarya/*_results/*
do
sed -i 's/FUN/MDAP/' $file
done


for file in ~/BE-OM-3_Bd-2-2_depl/CA.mr.49.17.15.0.02/funannotate_BE-OM-3.v1/*_results/*
do
sed -i 's/FUN/BDAP/' $file
done


for file in ~/Gv_SP/nextDenovo_polished_sgs-lgs_medaka/funannotate_Gv_SP.v1/*_results/*
do
sed -i 's/FUN/GDAP/' $file
done


for file in ~/assemblies/SK-39_Larssonia/megahit/funannotate_FI-SK-39.v1.1/*_results/*
do
sed -i 's/FUN/LOBT/' $file
done




cd /home/pascal/assemblies/comp_geno_Lobt

conda activate mamba-funannotate

export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/pasa-2.4.1

funannotate compare --cpus 12 --input ~/Gi_seq_second/funannotate_IL-YERU-16/annotate_results/Glugoides_intestinalis_IL-YERU-16.gbk \
~/assemblies/US/no_host/MAC_Ray81_metaspades/output/funannotate/US-FAR1-1/ \
~/PhD/data/trimmed_reads/spring_2017/SK-44_spr2017.megahit_asm/no_host/MAC_Ray81_metaspades/output/funannotate/SK-44/ \
~/assemblies/Mito_SP/masurca_ont_pb/CA.mr.99.17.15.0.02/contigs.SP.braker.dikarya/ \
~/assemblies/genebank/Paramicrosporidium_saccamoebae/ \
~/assemblies/genebank/Spraguea_lophii/ \
~/assemblies/genebank/Anncaliia_algerae/ \
~/assemblies/genebank/OC4/ \
~/assemblies/genebank/Encephalitozoon_romaleae/ \
~/assemblies/genebank/Nematocida_displodere/ \
~/assemblies/genebank/Nematocida_parisii \
~/assemblies/genebank/Nosema_ceranae/ \
~/assemblies/genebank/Rozella_allomycis/ \
~/assemblies/genebank/Vavraia_culicis_subsp._floridensis/ \
~/assemblies/genebank/Vittaforma_corneae/ \
~/assemblies/genebank/Edhazardia_aedis/ \
~/assemblies/genebank/Nosema_bombycis/ \
~/assemblies/genebank/Enterocytozoon_bieneusi/ \
~/assemblies/genebank/Hepatospora_eriocheir/ \
~/assemblies/genebank/Trachipleistophora_hominis/ \
~/assemblies/genebank/Pseudoloma_neurophilia/ \
~/assemblies/genebank/Hamiltosporidium_tvaerminnensis/ \
~/assemblies/genebank/Hamiltosporidium_magnivora/ \
~/assemblies/genebank/Tubulinosema_ratisbonensis/ \
~/assemblies/genebank/Dictyocoela_muelleri/ \
~/assemblies/genebank/Thelohania_contejeani/ \
~/assemblies/genebank/Cucumispora_dikerogammari/ \
~/BE-OM-3_Bd-2-2_depl/CA.mr.49.17.15.0.02/funannotate_BE-OM-3.v1/ \
~/Gv_SP/nextDenovo_polished_sgs-lgs_medaka/funannotate_Gv_SP.v2/ \
~/assemblies/SK-39_Larssonia/megahit/funannotate_FI-SK-39.v1.1/

#-------------------------------------------------------
#[Jun 13 01:43 PM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.12
#[Jun 13 01:43 PM]: Running 1.8.11
#[Jun 13 01:43 PM]: Now parsing 30 genomes
#[Jun 13 01:43 PM]: working on Glugoides intestinalis
#[Jun 13 01:43 PM]: working on Conglomerata sp.
#[Jun 13 01:43 PM]: working on Gurleya daphniae
#[Jun 13 01:43 PM]: working on Mitosporidium daphniae
#[Jun 13 01:43 PM]: working on Paramicrosporidium saccamoebae
#[Jun 13 01:43 PM]: working on Spraguea lophii 42_110
#[Jun 13 01:43 PM]: working on Anncaliia algerae PRA339
#[Jun 13 01:43 PM]: working on Ordospora colligata OC4
#[Jun 13 01:44 PM]: working on Encephalitozoon romaleae SJ-2008
#[Jun 13 01:44 PM]: working on Nematocida displodere
#[Jun 13 01:44 PM]: working on Nematocida parisii ERTm1
#[Jun 13 01:44 PM]: working on Nosema ceranae
#[Jun 13 01:44 PM]: working on Rozella allomycis CSF55
#[Jun 13 01:44 PM]: working on Vavraia culicis subsp. floridensis
#[Jun 13 01:44 PM]: working on Vittaforma corneae ATCC 50505
#[Jun 13 01:44 PM]: working on Edhazardia aedis USNM 41457
#[Jun 13 01:44 PM]: working on Nosema bombycis CQ1
#[Jun 13 01:44 PM]: working on Enterocytozoon bieneusi H348
#[Jun 13 01:44 PM]: working on Hepatospora eriocheir
#[Jun 13 01:44 PM]: working on Trachipleistophora hominis
#[Jun 13 01:44 PM]: working on Pseudoloma neurophilia
#[Jun 13 01:44 PM]: working on Hamiltosporidium tvaerminnensis
#[Jun 13 01:44 PM]: working on Hamiltosporidium magnivora
#[Jun 13 01:45 PM]: working on Tubulinosema ratisbonensis
#[Jun 13 01:45 PM]: working on Dictyocoela muelleri
#[Jun 13 01:45 PM]: working on Thelohania contejeani
#[Jun 13 01:45 PM]: working on Cucumispora dikerogammari
#[Jun 13 01:45 PM]: working on Binucleata daphniae
#[Jun 13 01:45 PM]: working on Gurleya daphniae
#[Jun 13 01:45 PM]: working on Conglomerata obtusa
#[Jun 13 01:45 PM]: No secondary metabolite annotations found
#[Jun 13 01:45 PM]: Summarizing PFAM domain results
#Fontconfig warning: ignoring UTF-8: not a valid region tag
#[Jun 13 01:45 PM]: Summarizing InterProScan results
#[Jun 13 01:45 PM]: Loading InterPro descriptions
#[Jun 13 01:45 PM]: Summarizing MEROPS protease results
#[Jun 13 01:46 PM]: found 16/87 MEROPS familes with stdev >= 1.000000
#[Jun 13 01:46 PM]: Summarizing CAZyme results
#[Jun 13 01:46 PM]: found 5/56 CAZy familes with stdev >= 1.000000
#[Jun 13 01:46 PM]: Summarizing COG results
#[Jun 13 01:46 PM]: Summarizing secreted protein results
#[Jun 13 01:46 PM]: Summarizing fungal transcription factors
#[Jun 13 01:46 PM]: Running GO enrichment for each genome
#  WARNING: skipping Hamiltosporidium_magnivora.txt as no GO terms
#  WARNING: skipping Edhazardia_aedis_USNM_41457.txt as no GO terms
#  WARNING: skipping Gurleya_daphniae.txt as no GO terms
#  WARNING: skipping Ordospora_colligata_OC4.txt as no GO terms
#  WARNING: skipping Nosema_bombycis_CQ1.txt as no GO terms
#  WARNING: skipping Nosema_ceranae.txt as no GO terms
#  WARNING: skipping Nematocida_displodere.txt as no GO terms
#  WARNING: skipping Trachipleistophora_hominis.txt as no GO terms
#  WARNING: skipping Vavraia_culicis_subsp._floridensis.txt as no GO terms
#  WARNING: skipping Rozella_allomycis_CSF55.txt as no GO terms
#  WARNING: skipping Paramicrosporidium_saccamoebae.txt as no GO terms
#  WARNING: skipping Hepatospora_eriocheir.txt as no GO terms
#  WARNING: skipping Vittaforma_corneae_ATCC_50505.txt as no GO terms
#  WARNING: skipping Hamiltosporidium_tvaerminnensis.txt as no GO terms
#  WARNING: skipping Anncaliia_algerae_PRA339.txt as no GO terms
#  WARNING: skipping Conglomerata_sp..txt as no GO terms
#  WARNING: skipping Thelohania_contejeani.txt as no GO terms
#  WARNING: skipping Spraguea_lophii_42_110.txt as no GO terms
#  WARNING: skipping Pseudoloma_neurophilia.txt as no GO terms
#  WARNING: skipping Tubulinosema_ratisbonensis.txt as no GO terms
#  WARNING: skipping Nematocida_parisii_ERTm1.txt as no GO terms
#  WARNING: skipping Enterocytozoon_bieneusi_H348.txt as no GO terms
#  WARNING: skipping Cucumispora_dikerogammari.txt as no GO terms
#  WARNING: skipping Encephalitozoon_romaleae_SJ-2008.txt as no GO terms
#  WARNING: skipping Dictyocoela_muelleri.txt as no GO terms
#/home/pascal/miniconda3/envs/mamba-funannotate/lib/python3.8/site-packages/funannotate/compare.py:804: FutureWarning: The default value of regex will change from True to False in a future version.
#  df.columns = df.columns.str.replace(r'^# ', '')
#[Jun 13 01:46 PM]: Running orthologous clustering tool, ProteinOrtho.  This may take awhile...
#[Jun 13 01:57 PM]: Compiling all annotations for each genome
#[Jun 13 01:58 PM]: Inferring phylogeny using iqtree
#[Jun 13 01:58 PM]: Found 4 single copy BUSCO orthologs, will use all to infer phylogeny
#[Jun 13 02:59 PM]: Compressing results to output file: funannotate_compare.tar.gz
#[Jun 13 03:01 PM]: Funannotate compare completed successfully!




