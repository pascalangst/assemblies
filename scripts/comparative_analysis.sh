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


## Glugoides hifi assembly

cd /home/pascal/assemblies/hifi/KZ-23-1



funannotate clean -i KZ-23-1.manual-merge.fa  -o KZ-23-1.clean.fa
funannotate sort -i KZ-23-1.clean.fa -o KZ-23-1.sort.fa
funannotate mask -i KZ-23-1.sort.fa -o KZ-23-1.mask.fa


# run funannotate with fungal gene hints
funannotate predict -i KZ-23-1.mask.fa -o funannotate_KZ-23-1 -s "Glugoides intestinalis" --isolate KZ-23-1 --augustus_species encephalitozoon_cuniculi_GB --busco_db microsporidia --busco_seed_species encephalitozoon_cuniculi_GB --cpus 12 \
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
#[Jul 12 03:08 PM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.13
#[Jul 12 03:08 PM]: Running funannotate v1.8.14
#[Jul 12 03:08 PM]: Skipping CodingQuarry as no --rna_bam passed
#[Jul 12 03:08 PM]: Parsed training data, run ab-initio gene predictors as follows:
#  Program      Training-Method
#  augustus     pretrained
#  genemark     selftraining
#  glimmerhmm   busco
#  snap         busco
#[Jul 12 03:08 PM]: Loading genome assembly and parsing soft-masked repetitive sequences
#[Jul 12 03:08 PM]: Genome loaded: 13 scaffolds; 5,125,750 bp; 5.70% repeats masked
#[Jul 12 03:08 PM]: Aligning transcript evidence to genome with minimap2
#[Jul 12 03:08 PM]: Found 15 alignments, wrote GFF3 and Augustus hints to file
#[Jul 12 03:08 PM]: Mapping 555,555 proteins to genome using diamond and exonerate
#[Jul 12 03:09 PM]: Found 96,433 preliminary alignments with diamond in 0:00:46 --> generated FASTA files for exonerate in 0:00:17
#     Progress: 96433 complete, 0 failed, 0 remaining
#[Jul 12 03:18 PM]: Exonerate finished in 0:09:00: found 22 alignments
#[Jul 12 03:18 PM]: Running GeneMark-ES on assembly
#[Jul 12 03:21 PM]: 2,487 predictions from GeneMark
#[Jul 12 03:21 PM]: Running BUSCO to find conserved gene models for training ab-initio predictors
#[Jul 12 03:24 PM]: 380 valid BUSCO predictions found, validating protein sequences
#[Jul 12 03:25 PM]: 378 BUSCO predictions validated
#[Jul 12 03:25 PM]: Running Augustus gene prediction using encephalitozoon_cuniculi_GB parameters
#     Progress: 16 complete, 0 failed, 0 remaining
#[Jul 12 03:26 PM]: 224 predictions from Augustus
#[Jul 12 03:26 PM]: Pulling out high quality Augustus predictions
#[Jul 12 03:26 PM]: Found 4 high quality predictions from Augustus (>90% exon evidence)
#[Jul 12 03:26 PM]: Running SNAP gene prediction, using training data: funannotate_KZ-23-1/predict_misc/busco.final.gff3
#[Jul 12 03:26 PM]: 2,256 predictions from SNAP
#[Jul 12 03:26 PM]: Running GlimmerHMM gene prediction, using training data: funannotate_KZ-23-1/predict_misc/busco.final.gff3
#[Jul 12 03:27 PM]: 2,736 predictions from GlimmerHMM
#[Jul 12 03:27 PM]: Summary of gene models passed to EVM (weights):
#  Source         Weight   Count
#  Augustus       1        220
#  Augustus HiQ   2        4
#  GeneMark       1        2487
#  GlimmerHMM     1        2736
#  snap           1        2256
#  Total          -        7703
#[Jul 12 03:27 PM]: EVM: partitioning input to ~ 35 genes per partition using min 1500 bp interval
#     Progress: 64 complete, 0 failed, 0 remaining
#[Jul 12 03:36 PM]: Converting to GFF3 and collecting all EVM results
#[Jul 12 03:36 PM]: 2,522 total gene models from EVM
#[Jul 12 03:36 PM]: Generating protein fasta files from 2,522 EVM models
#[Jul 12 03:36 PM]: now filtering out bad gene models (< 50 aa in length, transposable elements, etc).
#[Jul 12 03:36 PM]: Found 191 gene models to remove: 0 too short; 0 span gaps; 191 transposable elements
#[Jul 12 03:36 PM]: 2,331 gene models remaining
#[Jul 12 03:36 PM]: Predicting tRNAs
#[Jul 12 03:36 PM]: 51 tRNAscan models are valid (non-overlapping)
#[Jul 12 03:36 PM]: Generating GenBank tbl annotation file
#[Jul 12 03:36 PM]: Collecting final annotation files for 2,382 total gene models
#[Jul 12 03:36 PM]: Converting to final Genbank format
#[Jul 12 03:37 PM]: Funannotate predict is finished, output files are in the funannotate_KZ-23-1/predict_results folder
#[Jul 12 03:37 PM]: Your next step might be functional annotation, suggested commands:
#-------------------------------------------------------
#Run InterProScan (manual install):
#funannotate iprscan -i funannotate_KZ-23-1 -c 12
#
#Run antiSMASH (optional):
#funannotate remote -i funannotate_KZ-23-1 -m antismash -e youremail@server.edu
#
#Annotate Genome:
#funannotate annotate -i funannotate_KZ-23-1 --cpus 12 --sbt yourSBTfile.txt
#-------------------------------------------------------
#
#[Jul 12 03:37 PM]: Training parameters file saved: funannotate_KZ-23-1/predict_results/encephalitozoon_cuniculi_GB.parameters.json
#[Jul 12 03:37 PM]: Add species parameters to database:
#
#  funannotate species -s encephalitozoon_cuniculi_GB -a funannotate_KZ-23-1/predict_results/encephalitozoon_cuniculi_GB.parameters.json


mkdir funannotate_KZ-23-1/annotate_misc/
conda activate interproscan
python ~/Ht_methylation_project/analysis/iso-seq/long_read_protocol/iprscan-local.py -i funannotate_KZ-23-1/ -m local -c 12

conda activate mamba-funannotate-Ht
export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/pasa-2.4.1

funannotate annotate -i funannotate_KZ-23-1/ --cpus 12 --iprscan funannotate_KZ-23-1/annotate_misc/iprscan.xml --busco_db microsporidia
-------------------------------------------------------
#[Jul 12 04:16 PM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.13
#[Jul 12 04:16 PM]: Running 1.8.14
#[Jul 12 04:16 PM]: No NCBI SBT file given, will use default, however if you plan to submit to NCBI, create one and pass it here '--sbt'
#[Jul 12 04:16 PM]: Parsing input files
#[Jul 12 04:16 PM]: Existing tbl found: funannotate_KZ-23-1/predict_results/Glugoides_intestinalis_KZ-23-1.tbl
#[Jul 12 04:16 PM]: Adding Functional Annotation to Glugoides intestinalis, NCBI accession: None
#[Jul 12 04:16 PM]: Annotation consists of: 2,382 gene models
#[Jul 12 04:16 PM]: 2,331 protein records loaded
#[Jul 12 04:16 PM]: Running HMMer search of PFAM version 35.0
#[Jul 12 04:17 PM]: 1,897 annotations added
#[Jul 12 04:17 PM]: Running Diamond blastp search of UniProt DB version 2022_01
#[Jul 12 04:17 PM]: 41 valid gene/product annotations from 49 total
#[Jul 12 04:17 PM]: Running Eggnog-mapper
#[Jul 12 04:28 PM]: Parsing EggNog Annotations
#[Jul 12 04:28 PM]: EggNog version parsed as 2.1.9
#[Jul 12 04:28 PM]: 2,594 COG and EggNog annotations added
#[Jul 12 04:28 PM]: Combining UniProt/EggNog gene and product names using Gene2Product version 1.76
#[Jul 12 04:28 PM]: 619 gene name and product description annotations added
#[Jul 12 04:28 PM]: Running Diamond blastp search of MEROPS version 12.0
#[Jul 12 04:28 PM]: 48 annotations added
#[Jul 12 04:28 PM]: Annotating CAZYmes using HMMer search of dbCAN version 10.0
#[Jul 12 04:28 PM]: 6 annotations added
#[Jul 12 04:28 PM]: Annotating proteins with BUSCO microsporidia models
#[Jul 12 04:28 PM]: 515 annotations added
#[Jul 12 04:28 PM]: Predicting secreted and transmembrane proteins using Phobius
#     Progress: 2331 complete, 0 failed, 0 remaining
#[Jul 12 04:29 PM]: Predicting secreted proteins with SignalP
#[Jul 12 04:34 PM]: 100 secretome and 416 transmembane annotations added
#[Jul 12 04:34 PM]: Parsing InterProScan5 XML file
#[Jul 12 04:34 PM]: Found 0 duplicated annotations, adding 6,955 valid annotations
#[Jul 12 04:34 PM]: Converting to final Genbank format, good luck!
#[Jul 12 04:34 PM]: Creating AGP file and corresponding contigs file
#[Jul 12 04:34 PM]: Writing genome annotation table.
#[Jul 12 04:34 PM]: Funannotate annotate has completed successfully!
#
#        We need YOUR help to improve gene names/product descriptions:
#           0 gene/products names MUST be fixed, see funannotate_KZ-23-1/annotate_results/Gene2Products.must-fix.txt
#           6 gene/product names need to be curated, see funannotate_KZ-23-1/annotate_results/Gene2Products.need-curating.txt
#           11 gene/product names passed but are not in Database, see funannotate_KZ-23-1/annotate_results/Gene2Products.new-names-passed.txt
#
#        Please consider contributing a PR at https://github.com/nextgenusfs/gene2product
#
#-------------------------------------------------------





## Ordospora hifi assembly

cd ~/assemblies/hifi/GB-LK1-1/manual_merge/



funannotate clean -i GB-LK1-1.manual.merge.fa  -o GB-LK1-1.clean.fa
funannotate sort -i GB-LK1-1.clean.fa -o GB-LK1-1.sort.fa
funannotate mask -i GB-LK1-1.sort.fa -o GB-LK1-1.mask.fa


# run funannotate with fungal gene hints
funannotate predict -i GB-LK1-1.mask.fa -o funannotate_GB-LK1-1 -s "Ordospora colligata" --isolate GB-LK1-1 --augustus_species encephalitozoon_cuniculi_GB --busco_db microsporidia --busco_seed_species encephalitozoon_cuniculi_GB --cpus 12 \
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
#[Aug 07 01:54 PM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.13
#[Aug 07 01:54 PM]: Running funannotate v1.8.14
#
#[Aug 07 01:54 PM]: Skipping CodingQuarry as no --rna_bam passed
#[Aug 07 01:54 PM]: Parsed training data, run ab-initio gene predictors as follows:
#  Program      Training-Method
#  augustus     pretrained
#  genemark     selftraining
#  glimmerhmm   busco
#  snap         busco
#[Aug 07 01:54 PM]: Loading genome assembly and parsing soft-masked repetitive sequences
#[Aug 07 01:54 PM]: Genome loaded: 9 scaffolds; 3,181,619 bp; 11.93% repeats masked
#[Aug 07 01:54 PM]: Aligning transcript evidence to genome with minimap2
#[Aug 07 01:55 PM]: Found 1,867 alignments, wrote GFF3 and Augustus hints to file
#[Aug 07 01:55 PM]: Mapping 555,555 proteins to genome using diamond and exonerate
#[Aug 07 01:56 PM]: Found 99,116 preliminary alignments with diamond in 0:00:57 --> generated FASTA files for exonerate in 0:00:17
#     Progress: 99116 complete, 0 failed, 0 remaining
#[Aug 07 02:06 PM]: Exonerate finished in 0:10:04: found 335 alignments
#[Aug 07 02:06 PM]: Running GeneMark-ES on assembly
#[Aug 07 02:09 PM]: 1,393 predictions from GeneMark
#[Aug 07 02:09 PM]: Running BUSCO to find conserved gene models for training ab-initio predictors
#[Aug 07 02:12 PM]: 478 valid BUSCO predictions found, validating protein sequences
#[Aug 07 02:13 PM]: 477 BUSCO predictions validated
#[Aug 07 02:13 PM]: Running Augustus gene prediction using encephalitozoon_cuniculi_GB parameters
#     Progress: 9 complete, 0 failed, 0 remaining
#[Aug 07 02:13 PM]: 1,551 predictions from Augustus
#[Aug 07 02:13 PM]: Pulling out high quality Augustus predictions
#[Aug 07 02:13 PM]: Found 1,488 high quality predictions from Augustus (>90% exon evidence)
#[Aug 07 02:13 PM]: Running SNAP gene prediction, using training data: funannotate_GB-LK1-1/predict_misc/busco.final.gff3
#[Aug 07 02:14 PM]: 1,509 predictions from SNAP
#[Aug 07 02:14 PM]: Running GlimmerHMM gene prediction, using training data: funannotate_GB-LK1-1/predict_misc/busco.final.gff3
#[Aug 07 02:15 PM]: 1,825 predictions from GlimmerHMM
#[Aug 07 02:15 PM]: Summary of gene models passed to EVM (weights):
#  Source         Weight   Count
#  Augustus       1        66
#  Augustus HiQ   2        1497
#  GeneMark       1        1393
#  GlimmerHMM     1        1825
#  snap           1        1509
#  Total          -        6290
#[Aug 07 02:15 PM]: EVM: partitioning input to ~ 35 genes per partition using min 1500 bp interval
#     Progress: 23 complete, 0 failed, 0 remaining
#[Aug 07 02:30 PM]: Converting to GFF3 and collecting all EVM results
#[Aug 07 02:30 PM]: 1,774 total gene models from EVM
#[Aug 07 02:30 PM]: Generating protein fasta files from 1,774 EVM models
#[Aug 07 02:30 PM]: now filtering out bad gene models (< 50 aa in length, transposable elements, etc).
#[Aug 07 02:30 PM]: Found 6 gene models to remove: 0 too short; 0 span gaps; 6 transposable elements
#[Aug 07 02:30 PM]: 1,768 gene models remaining
#[Aug 07 02:30 PM]: Predicting tRNAs
#[Aug 07 02:31 PM]: 44 tRNAscan models are valid (non-overlapping)
#[Aug 07 02:31 PM]: Generating GenBank tbl annotation file
#[Aug 07 02:31 PM]: Collecting final annotation files for 1,812 total gene models
#[Aug 07 02:31 PM]: Converting to final Genbank format
#[Aug 07 02:31 PM]: Funannotate predict is finished, output files are in the funannotate_GB-LK1-1/predict_results folder
#[Aug 07 02:31 PM]: Your next step might be functional annotation, suggested commands:
#-------------------------------------------------------
#Run InterProScan (manual install):
#funannotate iprscan -i funannotate_GB-LK1-1 -c 12
#
#Run antiSMASH (optional):
#funannotate remote -i funannotate_GB-LK1-1 -m antismash -e youremail@server.edu
#
#Annotate Genome:
#funannotate annotate -i funannotate_GB-LK1-1 --cpus 12 --sbt yourSBTfile.txt
#-------------------------------------------------------
#
#[Aug 07 02:31 PM]: Training parameters file saved: funannotate_GB-LK1-1/predict_results/encephalitozoon_cuniculi_GB.parameters.json
#[Aug 07 02:31 PM]: Add species parameters to database:
#
#  funannotate species -s encephalitozoon_cuniculi_GB -a funannotate_GB-LK1-1/predict_results/encephalitozoon_cuniculi_GB.parameters.json

mkdir funannotate_GB-LK1-1/annotate_misc/
conda activate interproscan
python ~/Ht_methylation_project/analysis/iso-seq/long_read_protocol/iprscan-local.py -i funannotate_GB-LK1-1/ -m local -c 6

conda activate mamba-funannotate-Ht
export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/pasa-2.4.1

funannotate annotate -i funannotate_GB-LK1-1/ --cpus 12 --iprscan funannotate_GB-LK1-1/annotate_misc/iprscan.xml --busco_db microsporidia
#-------------------------------------------------------
#[Aug 07 05:15 PM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.13
#[Aug 07 05:15 PM]: Running 1.8.14
#[Aug 07 05:15 PM]: No NCBI SBT file given, will use default, however if you plan to submit to NCBI, create one and pass it here '--sbt'
#[Aug 07 05:15 PM]: Parsing input files
#[Aug 07 05:15 PM]: Existing tbl found: funannotate_GB-LK1-1/predict_results/Ordospora_colligata_GB-LK1-1.tbl
#[Aug 07 05:15 PM]: Adding Functional Annotation to Ordospora colligata, NCBI accession: None
#[Aug 07 05:15 PM]: Annotation consists of: 1,812 gene models
#[Aug 07 05:15 PM]: 1,768 protein records loaded
#[Aug 07 05:15 PM]: Running HMMer search of PFAM version 35.0
#[Aug 07 05:16 PM]: 1,754 annotations added
#[Aug 07 05:16 PM]: Running Diamond blastp search of UniProt DB version 2022_01
#[Aug 07 05:17 PM]: 160 valid gene/product annotations from 223 total
#[Aug 07 05:17 PM]: Running Eggnog-mapper
#[Aug 07 05:30 PM]: Parsing EggNog Annotations
#[Aug 07 05:30 PM]: EggNog version parsed as 2.1.9
#[Aug 07 05:30 PM]: 2,546 COG and EggNog annotations added
#[Aug 07 05:30 PM]: Combining UniProt/EggNog gene and product names using Gene2Product version 1.76
#[Aug 07 05:30 PM]: 650 gene name and product description annotations added
#[Aug 07 05:30 PM]: Running Diamond blastp search of MEROPS version 12.0
#[Aug 07 05:30 PM]: 42 annotations added
#[Aug 07 05:30 PM]: Annotating CAZYmes using HMMer search of dbCAN version 10.0
#[Aug 07 05:30 PM]: 7 annotations added
#[Aug 07 05:30 PM]: Annotating proteins with BUSCO microsporidia models
#[Aug 07 05:30 PM]: 483 annotations added
#[Aug 07 05:30 PM]: Predicting secreted and transmembrane proteins using Phobius
#     Progress: 1768 complete, 0 failed, 0 remaining
#[Aug 07 05:30 PM]: Predicting secreted proteins with SignalP
#[Aug 07 05:36 PM]: 57 secretome and 328 transmembane annotations added
#[Aug 07 05:36 PM]: Parsing InterProScan5 XML file
#[Aug 07 05:36 PM]: Found 0 duplicated annotations, adding 6,637 valid annotations
#[Aug 07 05:36 PM]: Converting to final Genbank format, good luck!
#[Aug 07 05:36 PM]: Creating AGP file and corresponding contigs file
#[Aug 07 05:36 PM]: Writing genome annotation table.
#[Aug 07 05:36 PM]: Funannotate annotate has completed successfully!
#
#        We need YOUR help to improve gene names/product descriptions:
#           0 gene/products names MUST be fixed, see funannotate_GB-LK1-1/annotate_results/Gene2Products.must-fix.txt
#           4 gene/product names need to be curated, see funannotate_GB-LK1-1/annotate_results/Gene2Products.need-curating.txt
#           15 gene/product names passed but are not in Database, see funannotate_GB-LK1-1/annotate_results/Gene2Products.new-names-passed.txt
#
#        Please consider contributing a PR at https://github.com/nextgenusfs/gene2product

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


for file in ~/assemblies/hifi/GB-LK1-1/manual_merge/funannotate_GB-LK1-1/*_results/*
do
sed -i 's/FUN/OCOL/' $file
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






for file in ~/Gv_SP/nextDenovo_polished_sgs-lgs_medaka/funannotate_Gv_SP.v2/*_results/*
do
sed -i 's/FUN/GVAV/' $file
done



cd /home/pascal/assemblies/comp_geno_PhD_Vnec_GiHiFi

conda activate mamba-funannotate

export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/pasa-2.4.1

funannotate compare --cpus 18 --input ~/assemblies/hifi/KZ-23-1/funannotate_KZ-23-1/ \
~/assemblies/Mito_SP/masurca_ont_pb/CA.mr.99.17.15.0.02/contigs.SP.braker.dikarya/ \
~/assemblies/genebank/Paramicrosporidium_saccamoebae/ \
~/assemblies/genebank/Spraguea_lophii/ \
~/assemblies/genebank/Anncaliia_algerae/ \
~/assemblies/genebank/Encephalitozoon_intestinalis/ \
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
~/assemblies/genebank/Hamiltosporidium_magnivora/ \
~/assemblies/genebank/Tubulinosema_ratisbonensis/ \
~/assemblies/genebank/Dictyocoela_muelleri/ \
~/assemblies/genebank/Thelohania_contejeani/ \
~/assemblies/genebank/Cucumispora_dikerogammari/ \
~/BE-OM-3_Bd-2-2_depl/CA.mr.49.17.15.0.02/funannotate_BE-OM-3.v1/ \
~/Gv_SP/nextDenovo_polished_sgs-lgs_medaka/funannotate_Gv_SP.v2/ \
~/assemblies/SK-39_Larssonia/megahit/funannotate_FI-SK-39.v1.1/ \
~/assemblies/hifi/GB-LK1-1/manual_merge/funannotate_GB-LK1-1/ \
~/Ht_methylation_project/analysis/iso-seq/long_read_protocol/FI-OER-3-3.v2.6/ \
~/assemblies/genebank/Vairimorpha_necatrix/
#-------------------------------------------------------
#[12/20/24 14:57:43]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.12
#[12/20/24 14:57:43]: Running 1.8.11
#[12/20/24 14:57:49]: Now parsing 29 genomes
#[12/20/24 14:57:49]: working on Glugoides intestinalis
#[12/20/24 14:57:52]: working on Mitosporidium daphniae
#[12/20/24 14:57:56]: working on Paramicrosporidium saccamoebae
#[12/20/24 14:58:00]: working on Spraguea lophii 42_110
#[12/20/24 14:58:03]: working on Anncaliia algerae PRA339
#[12/20/24 14:58:06]: working on Encephalitozoon intestinalis
#[12/20/24 14:58:08]: working on Nematocida displodere
#[12/20/24 14:58:10]: working on Nematocida parisii ERTm1
#[12/20/24 14:58:13]: working on Nosema ceranae
#[12/20/24 14:58:17]: working on Rozella allomycis CSF55
#[12/20/24 14:58:26]: working on Vavraia culicis subsp. floridensis
#[12/20/24 14:58:28]: working on Vittaforma corneae ATCC 50505
#[12/20/24 14:58:31]: working on Edhazardia aedis USNM 41457
#[12/20/24 14:58:38]: working on Nosema bombycis CQ1
#[12/20/24 14:58:44]: working on Enterocytozoon bieneusi H348
#[12/20/24 14:58:48]: working on Hepatospora eriocheir
#[12/20/24 14:58:51]: working on Trachipleistophora hominis
#[12/20/24 14:58:54]: working on Pseudoloma neurophilia
#[12/20/24 14:58:59]: working on Hamiltosporidium magnivora
#[12/20/24 14:59:05]: working on Tubulinosema ratisbonensis
#[12/20/24 14:59:11]: working on Dictyocoela muelleri
#[12/20/24 14:59:23]: working on Thelohania contejeani
#[12/20/24 14:59:28]: working on Cucumispora dikerogammari
#[12/20/24 14:59:37]: working on Binucleata daphniae
#[12/20/24 14:59:40]: working on Gurleya daphniae
#[12/20/24 14:59:45]: working on Conglomerata obtusa
#[12/20/24 14:59:51]: working on Ordospora colligata
#[12/20/24 14:59:54]: working on Hamiltosporidium tvaerminnensis
#[12/20/24 14:59:58]: working on Vairimorpha necatrix
#[12/20/24 15:00:01]: No secondary metabolite annotations found
#[12/20/24 15:00:01]: Summarizing PFAM domain results
#[12/20/24 15:00:12]: Summarizing InterProScan results
#[12/20/24 15:00:12]: Loading InterPro descriptions
#[12/20/24 15:00:12]: Summarizing MEROPS protease results
#[12/20/24 15:00:13]: found 17/87 MEROPS familes with stdev >= 1.000000
#[12/20/24 15:00:14]: Summarizing CAZyme results
#[12/20/24 15:00:14]: found 4/56 CAZy familes with stdev >= 1.000000
#[12/20/24 15:00:15]: Summarizing COG results
#12/20/24 15:00:20]: Summarizing secreted protein results
#[12/20/24 15:00:21]: Summarizing fungal transcription factors
#[12/20/24 15:00:23]: Running GO enrichment for each genome
#  WARNING: skipping Hamiltosporidium_magnivora.txt as no GO terms
#  WARNING: skipping Edhazardia_aedis_USNM_41457.txt as no GO terms
#  WARNING: skipping Encephalitozoon_intestinalis.txt as no GO terms
#  WARNING: skipping Nosema_bombycis_CQ1.txt as no GO terms
#  WARNING: skipping Nosema_ceranae.txt as no GO terms
#  WARNING: skipping Nematocida_displodere.txt as no GO terms
#  WARNING: skipping Trachipleistophora_hominis.txt as no GO terms
#  WARNING: skipping Vavraia_culicis_subsp._floridensis.txt as no GO terms
#  WARNING: skipping Rozella_allomycis_CSF55.txt as no GO terms
#  WARNING: skipping Paramicrosporidium_saccamoebae.txt as no GO terms
#  WARNING: skipping Hepatospora_eriocheir.txt as no GO terms
#  WARNING: skipping Vittaforma_corneae_ATCC_50505.txt as no GO terms
#  WARNING: skipping Anncaliia_algerae_PRA339.txt as no GO terms
#  WARNING: skipping Vairimorpha_necatrix.txt as no GO terms
#  WARNING: skipping Thelohania_contejeani.txt as no GO terms
#  WARNING: skipping Spraguea_lophii_42_110.txt as no GO terms
#  WARNING: skipping Pseudoloma_neurophilia.txt as no GO terms
#  WARNING: skipping Tubulinosema_ratisbonensis.txt as no GO terms
#  WARNING: skipping Nematocida_parisii_ERTm1.txt as no GO terms
#  WARNING: skipping Enterocytozoon_bieneusi_H348.txt as no GO terms
#  WARNING: skipping Cucumispora_dikerogammari.txt as no GO terms
#  WARNING: skipping Dictyocoela_muelleri.txt as no GO terms
#/home/pascal/miniconda3/envs/mamba-funannotate/lib/python3.8/site-packages/funannotate/compare.py:804: FutureWarning: The default value of regex will change from True to False in a future #version.
#  df.columns = df.columns.str.replace(r'^# ', '')
#[12/20/24 15:00:59]: Running orthologous clustering tool, ProteinOrtho.  This may take awhile...
#[12/20/24 15:10:33]: Compiling all annotations for each genome
#[12/20/24 15:11:03]: Inferring phylogeny using iqtree
#[12/20/24 15:11:03]: Found 6 single copy BUSCO orthologs, will use all to infer phylogeny
#[12/20/24 18:12:22]: Compressing results to output file: funannotate_compare.tar.gz
#[12/20/24 18:14:44]: Funannotate compare completed successfully!

# do phylogeny based on the identified SCOs
proteinortho_grab_proteins.pl -exact 'H312_01630-T1,BDAP_002774-T1,COBT_003363-T1,CDIK_1150-T1,DMUE_2113-T1,EDEG_01226-T1,GPK93_04g05510-T1,EBI_22833-T1,GINT2_000391-T1,GVAV_001470-T1,CWI36_3044p0010-T1,FUN_000534-T1,HERIO_1446-T1,MDAP_000294-T1,NEDG_00316-T1,NEPG_00291-T1,NBO_554g0006-T1,G9O61_00g011750-T1,OCOL_001001-T1,PSACC_03723-T1,M153_1100042559-T1,ROZALSC1DRAFT_27770-T1,SLOPH_25-T1,TCON_1062-T1,THOM_1010-T1,TUBRATIS_000190-T1,VNE69_01343-T1,VCUG_00901-T1,VICG_00544-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 1.fasta
proteinortho_grab_proteins.pl -exact 'H312_00699-T1,BDAP_001326-T1,COBT_001972-T1,CDIK_1973-T1,DMUE_2428-T1,EDEG_03183-T1,GPK93_02g01900-T1,EBI_21688-T1,GINT2_000016-T1,GVAV_001045-T1,CWI36_0678p0010-T1,FUN_003080-T1,HERIO_1499-T1,MDAP_001592-T1,NEDG_00370-T1,NEPG_00197-T1,NBO_4g0069-T1,G9O61_00g018800-T1,OCOL_001361-T1,PSACC_00668-T1,M153_13560001237-T1,ROZALSC1DRAFT_29250-T1,SLOPH_7-T1,TCON_0101-T1,THOM_2344-T1,TUBRATIS_22680-T1,VNE69_06064-T1,VCUG_01595-T1,VICG_00714-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 2.fasta
proteinortho_grab_proteins.pl -exact 'H312_02923-T1,BDAP_000450-T1,COBT_002239-T1,CDIK_1788-T1,DMUE_5550-T1,EDEG_02940-T1,GPK93_11g21110-T1,EBI_24670-T1,GINT2_001882-T1,GVAV_002207-T1,CWI36_0599p0050-T1,FUN_3534-T1,HERIO_924-T1,MDAP_002728-T1,NEDG_00805-T1,NEPG_01108-T1,NBO_10g0045-T1,G9O61_00g002220-T1,OCOL_001245-T1,PSACC_01566-T1,M153_6500024225-T1,ROZALSC1DRAFT_28270-T1,SLOPH_1999-T1,TCON_2475-T1,THOM_3138-T1,TUBRATIS_001240-T1,VNE69_01127-T1,VCUG_01638-T1,VICG_01056-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 3.fasta
proteinortho_grab_proteins.pl -exact 'H312_00354-T1,BDAP_001113-T1,COBT_001757-T1,CDIK_0650-T1,DMUE_4803-T1,EDEG_02631-T1,GPK93_10g17830-T1,EBI_27246-T1,GINT2_001092-T1,GVAV_000067-T1,CWI36_0540p0020-T1,FUN_001200-T1,HERIO_2105-T1,MDAP_002690-T1,NEDG_02160-T1,NEPG_02164-T1,NBO_69g0011-T1,G9O61_00g004770-T1,OCOL_000184-T1,PSACC_03103-T1,M153_3860006633-T1,ROZALSC1DRAFT_26648-T1,SLOPH_1230-T1,TCON_0396-T1,THOM_3117-T1,TUBRATIS_26120-T1,VNE69_07248-T1,VCUG_00398-T1,VICG_01562-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 4.fasta
proteinortho_grab_proteins.pl -exact 'H312_01305-T1,BDAP_001359-T1,COBT_000603-T1,CDIK_2702-T1,DMUE_1538-T1,EDEG_01201-T1,GPK93_04g06340-T1,EBI_27095-T1,GINT2_001633-T1,GVAV_002143-T1,CWI36_1861p0020-T1,FUN_002515-T1,HERIO_1181-T1,MDAP_001741-T1,NEDG_00621-T1,NEPG_00865-T1,NBO_32g0045-T1,G9O61_00g009880-T1,OCOL_001069-T1,PSACC_00724-T1,M153_1570003290-T1,ROZALSC1DRAFT_28387-T1,SLOPH_345-T1,TCON_0250-T1,THOM_1775-T1,TUBRATIS_003730-T1,VNE69_12125-T1,VCUG_00672-T1,VICG_00792-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 5.fasta
proteinortho_grab_proteins.pl -exact 'H312_00789-T1,BDAP_002877-T1,COBT_002788-T1,CDIK_3368-T1,DMUE_0402-T1,EDEG_01293-T1,GPK93_05g08030-T1,EBI_27165-T1,GINT2_001595-T1,GVAV_001226-T1,CWI36_0280p0020-T1,FUN_002223-T1,HERIO_1310-T1,MDAP_002521-T1,NEDG_00128-T1,NEPG_00519-T1,NBO_72g0004-T1,G9O61_00g016100-T1,OCOL_000621-T1,PSACC_01191-T1,M153_3560004869-T1,ROZALSC1DRAFT_27522-T1,SLOPH_48-T1,TCON_2510-T1,THOM_2953-T1,TUBRATIS_13470-T1,VNE69_05175-T1,VCUG_02140-T1,VICG_01910-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 6.fasta
proteinortho_grab_proteins.pl -exact 'H312_03235-T1,BDAP_001452-T1,COBT_000670-T1,CDIK_0066-T1,DMUE_2713-T1,EDEG_01612-T1,GPK93_10g17700-T1,EBI_27253-T1,GINT2_001103-T1,GVAV_001478-T1,CWI36_0560p0040-T1,FUN_001501-T1,HERIO_1929-T1,MDAP_001577-T1,NEDG_01648-T1,NEPG_02004-T1,NBO_507g0006-T1,G9O61_00g004630-T1,OCOL_000197-T1,PSACC_02324-T1,M153_1364000159-T1,ROZALSC1DRAFT_29961-T1,SLOPH_2456-T1,TCON_1274-T1,THOM_1444-T1,TUBRATIS_16190-T1,VNE69_07264-T1,VCUG_00041-T1,VICG_01546-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 7.fasta
proteinortho_grab_proteins.pl -exact 'H312_02747-T1,BDAP_002071-T1,COBT_001244-T1,CDIK_2281-T1,DMUE_2984-T1,EDEG_02789-T1,GPK93_07g11480-T1,EBI_24199-T1,GINT2_001989-T1,GVAV_003518-T1,CWI36_1856p0020-T1,FUN_001997-T1,HERIO_768-T1,MDAP_000377-T1,NEDG_00049-T1,NEPG_00578-T1,NBO_34g0010-T1,G9O61_00g006930-T1,OCOL_000815-T1,PSACC_02449-T1,M153_7000003463-T1,ROZALSC1DRAFT_19879-T1,SLOPH_1358-T1,TCON_0485-T1,THOM_0830-T1,TUBRATIS_006320-T1,VNE69_03308-T1,VCUG_02629-T1,VICG_00363-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 8.fasta
proteinortho_grab_proteins.pl -exact 'H312_00419-T1,BDAP_000210-T1,COBT_001983-T1,CDIK_0550-T1,DMUE_1251-T1,EDEG_02786-T1,GPK93_10g18270-T1,EBI_22821-T1,GINT2_000077-T1,GVAV_000779-T1,CWI36_0568p0030-T1,FUN_000329-T1,HERIO_2138-T1,MDAP_000799-T1,NEDG_01364-T1,NEPG_01782-T1,NBO_375g0006-T1,G9O61_00g005330-T1,OCOL_000145-T1,PSACC_00101-T1,M153_13690002889-T1,ROZALSC1DRAFT_27300-T1,SLOPH_2461-T1,TCON_0158-T1,THOM_2146-T1,TUBRATIS_002770-T1,VNE69_07186-T1,VCUG_00546-T1,VICG_00073-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 9.fasta
proteinortho_grab_proteins.pl -exact 'H312_03238-T1,BDAP_001643-T1,COBT_000402-T1,CDIK_0806-T1,DMUE_0650-T1,EDEG_01378-T1,GPK93_05g08450-T1,EBI_23231-T1,GINT2_001648-T1,GVAV_003468-T1,CWI36_0114p0040-T1,FUN_001927-T1,HERIO_27-T1,MDAP_000989-T1,NEDG_00898-T1,NEPG_02477-T1,NBO_609gi001-T1,G9O61_00g000360-T1,OCOL_000032-T1,PSACC_00251-T1,M153_4127000461-T1,ROZALSC1DRAFT_29142-T1,SLOPH_2294-T1,TCON_2685-T1,THOM_3095-T1,TUBRATIS_23980-T1,VNE69_08047-T1,VCUG_01731-T1,VICG_00597-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 10.fasta
proteinortho_grab_proteins.pl -exact 'H312_02256-T1,BDAP_002336-T1,COBT_001645-T1,CDIK_2229-T1,DMUE_4804-T1,EDEG_00716-T1,GPK93_05g07860-T1,EBI_24873-T1,GINT2_000825-T1,GVAV_000770-T1,CWI36_0215p0030-T1,FUN_002921-T1,HERIO_521-T1,MDAP_000206-T1,NEDG_01458-T1,NEPG_01244-T1,NBO_53g0029-T1,G9O61_00g009230-T1,OCOL_000635-T1,PSACC_02963-T1,M153_2200020223-T1,ROZALSC1DRAFT_31056-T1,SLOPH_310-T1,TCON_0317-T1,THOM_1380-T1,TUBRATIS_008340-T1,VNE69_05117-T1,VCUG_00209-T1,VICG_01603-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 11.fasta
proteinortho_grab_proteins.pl -exact 'H312_01166-T1,BDAP_000524-T1,COBT_001446-T1,CDIK_0267-T1,DMUE_2627-T1,EDEG_01657-T1,GPK93_04g05720-T1,EBI_22378-T1,GINT2_000379-T1,GVAV_000581-T1,CWI36_1222p0020-T1,FUN_000501-T1,HERIO_61-T1,MDAP_000302-T1,NEDG_01152-T1,NEPG_01329-T1,NBO_64g0025-T1,G9O61_00g012100-T1,OCOL_001015-T1,PSACC_03732-T1,M153_2500035249-T1,ROZALSC1DRAFT_27593-T1,SLOPH_1701-T1,TCON_2572-T1,THOM_0227-T1,TUBRATIS_003500-T1,VNE69_01382-T1,VCUG_00940-T1,VICG_00509-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 12.fasta
proteinortho_grab_proteins.pl -exact 'H312_00444-T1,BDAP_002682-T1,COBT_000940-T1,CDIK_1276-T1,DMUE_2053-T1,EDEG_03573-T1,GPK93_10g18350-T1,EBI_26123-T1,GINT2_000052-T1,GVAV_000578-T1,CWI36_0328p0040-T1,FUN_000275-T1,HERIO_1880-T1,MDAP_000158-T1,NEDG_01329-T1,NEPG_01739-T1,NBO_84g0001-T1,G9O61_00g005280-T1,OCOL_000138-T1,PSACC_00497-T1,M153_1890003809-T1,ROZALSC1DRAFT_29469-T1,SLOPH_2478-T1,TCON_0145-T1,THOM_2264-T1,TUBRATIS_12430-T1,VNE69_07190-T1,VCUG_01152-T1,VICG_00081-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 13.fasta
proteinortho_grab_proteins.pl -exact 'H312_00416-T1,BDAP_001806-T1,COBT_001189-T1,CDIK_0516-T1,DMUE_4879-T1,EDEG_00441-T1,GPK93_08g13090-T1,EBI_24859-T1,GINT2_001945-T1,GVAV_001327-T1,CWI36_2446p0020-T1,FUN_001473-T1,HERIO_382-T1,MDAP_001337-T1,NEDG_01757-T1,NEPG_01907-T1,NBO_59g0003-T1,G9O61_00g017520-T1,OCOL_001652-T1,PSACC_03357-T1,M153_27110001801-T1,ROZALSC1DRAFT_6479-T1,SLOPH_2200-T1,TCON_0245-T1,THOM_1918-T1,TUBRATIS_16600-T1,VNE69_02031-T1,VCUG_02635-T1,VICG_01804-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 14.fasta
proteinortho_grab_proteins.pl -exact 'H312_02108-T1,BDAP_001437-T1,COBT_000551-T1,CDIK_0540-T1,DMUE_2398-T1,EDEG_02744-T1,GPK93_08g14520-T1,EBI_21958-T1,GINT2_002329-T1,GVAV_000137-T1,CWI36_0156p0060-T1,FUN_000113-T1,HERIO_789-T1,MDAP_002125-T1,NEDG_01888-T1,NEPG_01597-T1,NBO_12g0011-T1,G9O61_00g022470-T1,OCOL_001771-T1,PSACC_02817-T1,M153_7060003767-T1,ROZALSC1DRAFT_30734-T1,SLOPH_2047-T1,TCON_0242-T1,THOM_1367-T1,TUBRATIS_16960-T1,VNE69_10031-T1,VCUG_00195-T1,VICG_00889-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_intestinalis.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vairimorpha_necatrix.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 15.fasta




#sorted.species.list
#>Binucleata daphniae
#>Cucumispora dikerogammari
#>Conglomerata obtusa
#>Hamiltosporidium magnivora
#>Dictyocoela muelleri
#>Enterocytozoon bieneusi
#>Edhazardia aedis
#>Hamiltosporidium tvaerminnensis
#>Nosema ceranae
#>Glugoides intestinalis
#>Encephalitozoon intestinalis
#>Gurleya vavrai
#>Anncaliia algerae
#>Hepatospora eriocheir
#>Pseudoloma neurophilia
#>Mitosporidium daphniae
#>Nosema bombycis
#>Nematocida displodere
#>Nematocida parisii
#>Ordospora colligata
#>Paramicrosporidium saccamoebae
#>Rozella allomycis
#>Spraguea lophii
#>Astathelohania contejeani
#>Trachipleistophora hominis
#>Tubulinosema ratisbonensis
#>Vavraia culicis subsp. floridensis
#>Vittaforma corneae
#>Vairimorpha necatrix


# sort sequences by header and replace header with species name
for fasta in *.fasta
do
paste -d'\n' sorted.species.list <(awk 'NR%2==0' <(bioawk -c fastx '{print}' $fasta | sort -k1,1V | awk '{print ">"$1;print $2}')) > $fasta.sorted.fasta
done

# infer gene trees from SCOs (align, trim, make tree)
for fasta in *sorted.fasta
do
mafft --anysymbol --quiet $fasta > $fasta.mafft.fasta
sed -i "" 's/ /_/g' $fasta.mafft.fasta
trimal -in $fasta.mafft.fasta -out $fasta.mafft.trimal.phylip -automated1 -phylip
trimal -in $fasta.mafft.fasta -out $fasta.mafft.trimal.fasta -automated1 -fasta
iqtree2 -s $fasta.mafft.trimal.phylip -nt AUTO -ntmax 12 -seed 12345 -bb 1000 -o "Rozella_allomycis"
done

# ASTRAL "weigthed" phylogeny based on gene trees
cat *fasta.sorted.fasta.mafft.trimal.phylip.treefile > all.fasta.sorted.fasta.mafft.trimal.phylip.treefile
~/bioinformatics/ASTER-MacOS/bin/wastral -r 16 -s 16 -x 100 -n 0 --root "Rozella_allomycis" -t 12 -o all.fasta.sorted.fasta.mafft.trimal.phylip.astral all.fasta.sorted.fasta.mafft.trimal.phylip.treefile 2>all.fasta.sorted.fasta.mafft.trimal.phylip.log

# concatenate all aligned and trimmed SCOs
#seqkit concat *sorted.fasta | sed -e 's/\|.*$//' > sorted.concat.fasta
cut -f1 *.trimal.fasta | awk -v RS=">" -v FS="\n" -v OFS="\n" '{for(i=2; i<=NF; i++) {seq[$1] = seq[$1]$i}}; END {for(id in seq){print ">"id, seq[id]}}' > sorted.mafft.trimal.concat.fasta

# infer phylogeny
iqtree2 -s sorted.mafft.trimal.concat.fasta -nt AUTO -ntmax 12 -seed 12345 -bb 1000 -o "Rozella_allomycis"







### Focusing on the studied species only
cd /home/pascal/assemblies/comp_geno_PhD_reduced

conda activate mamba-funannotate

export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/pasa-2.4.1

funannotate compare --cpus 18 --input ~/assemblies/hifi/KZ-23-1/funannotate_KZ-23-1/ \
~/assemblies/Mito_SP/masurca_ont_pb/CA.mr.99.17.15.0.02/contigs.SP.braker.dikarya/ \
~/assemblies/genebank/Encephalitozoon_intestinalis/ \
~/BE-OM-3_Bd-2-2_depl/CA.mr.49.17.15.0.02/funannotate_BE-OM-3.v1/ \
~/Gv_SP/nextDenovo_polished_sgs-lgs_medaka/funannotate_Gv_SP.v2/ \
~/assemblies/SK-39_Larssonia/megahit/funannotate_FI-SK-39.v1.1/ \
~/assemblies/hifi/GB-LK1-1/manual_merge/funannotate_GB-LK1-1/ \
~/Ht_methylation_project/analysis/iso-seq/long_read_protocol/FI-OER-3-3.v2.6/ \
~/assemblies/genebank/Vairimorpha_necatrix/
#-------------------------------------------------------
#[Dec 20 06:14 PM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.12
#[Dec 20 06:14 PM]: Running 1.8.11
#[Dec 20 06:14 PM]: Now parsing 9 genomes
#[Dec 20 06:14 PM]: working on Glugoides intestinalis
#[Dec 20 06:14 PM]: working on Mitosporidium daphniae
#[Dec 20 06:14 PM]: working on Encephalitozoon intestinalis
#[Dec 20 06:14 PM]: working on Binucleata daphniae
#[Dec 20 06:14 PM]: working on Gurleya daphniae
#[Dec 20 06:15 PM]: working on Conglomerata obtusa
#[Dec 20 06:15 PM]: working on Ordospora colligata
#[Dec 20 06:15 PM]: working on Hamiltosporidium tvaerminnensis
#[Dec 20 06:15 PM]: working on Vairimorpha necatrix
#[Dec 20 06:15 PM]: No secondary metabolite annotations found
#[Dec 20 06:15 PM]: Summarizing PFAM domain results
#Fontconfig warning: ignoring UTF-8: not a valid region tag
#[Dec 20 06:15 PM]: Summarizing InterProScan results
#[Dec 20 06:15 PM]: Loading InterPro descriptions
#[Dec 20 06:15 PM]: Summarizing MEROPS protease results
#[Dec 20 06:15 PM]: found 11/52 MEROPS familes with stdev >= 1.000000
#[Dec 20 06:15 PM]: Summarizing CAZyme results
#[Dec 20 06:15 PM]: found 19 CAZy familes
#[Dec 20 06:15 PM]: Summarizing COG results
#[Dec 20 06:15 PM]: Summarizing secreted protein results
#[Dec 20 06:15 PM]: Summarizing fungal transcription factors
#[Dec 20 06:15 PM]: Running GO enrichment for each genome
#  WARNING: skipping Encephalitozoon_intestinalis.txt as no GO terms
#  WARNING: skipping Vairimorpha_necatrix.txt as no GO terms
#/home/pascal/miniconda3/envs/mamba-funannotate/lib/python3.8/site-packages/funannotate/compare.py:804: FutureWarning: The default value of regex will change from True to False in a future version.
#  df.columns = df.columns.str.replace(r'^# ', '')
#[Dec 20 06:16 PM]: Running orthologous clustering tool, ProteinOrtho.  This may take awhile...
#[Dec 20 06:17 PM]: Compiling all annotations for each genome
#[Dec 20 06:17 PM]: Inferring phylogeny using iqtree
#[Dec 20 06:17 PM]: Found 73 single copy BUSCO orthologs, will use all to infer phylogeny
#[Dec 20 08:21 PM]: Compressing results to output file: funannotate_compare.tar.gz
#[Dec 20 08:22 PM]: Funannotate compare completed successfully!

# do phylogeny based on the identified SCOs
proteinortho_grab_proteins.pl -exact 'BDAP_000413-T1,COBT_001743-T1,GPK93_11g20160-T1,GINT2_000264-T1,GVAV_002274-T1,FUN_000773-T1,MDAP_001513-T1,OCOL_001168-T1,VNE69_01196-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 1.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002683-T1,COBT_000941-T1,GPK93_10g18340-T1,GINT2_000073-T1,GVAV_000579-T1,FUN_000261-T1,MDAP_002268-T1,OCOL_000139-T1,VNE69_07192-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 2.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000317-T1,COBT_001449-T1,GPK93_11g20860-T1,GINT2_000984-T1,GVAV_001294-T1,FUN_000665-T1,MDAP_001208-T1,OCOL_001221-T1,VNE69_01260-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 3.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001331-T1,COBT_001454-T1,GPK93_04g06680-T1,GINT2_002360-T1,GVAV_002162-T1,FUN_002466-T1,MDAP_002068-T1,OCOL_001096-T1,VNE69_11152-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 4.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001069-T1,COBT_002408-T1,GPK93_09g17160-T1,GINT2_001453-T1,GVAV_000165-T1,FUN_001052-T1,MDAP_001747-T1,OCOL_000446-T1,VNE69_03234-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 5.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001549-T1,COBT_002425-T1,GPK93_05g07450-T1,GINT2_000552-T1,GVAV_001032-T1,FUN_003086-T1,MDAP_000677-T1,OCOL_000673-T1,VNE69_05053-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 6.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002819-T1,COBT_000242-T1,GPK93_03g04660-T1,GINT2_000484-T1,GVAV_001502-T1,FUN_000503-T1,MDAP_000513-T1,OCOL_000252-T1,VNE69_01390-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 7.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001639-T1,COBT_002683-T1,GPK93_06g09610-T1,GINT2_000334-T1,GVAV_003455-T1,FUN_002694-T1,MDAP_001529-T1,OCOL_001408-T1,VNE69_01033-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 8.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002048-T1,COBT_000494-T1,GPK93_07g10870-T1,GINT2_002003-T1,GVAV_000374-T1,FUN_001907-T1,MDAP_002223-T1,OCOL_001595-T1,VNE69_03334-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 9.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002896-T1,COBT_001430-T1,GPK93_04g05800-T1,GINT2_002260-T1,GVAV_000509-T1,FUN_001161-T1,MDAP_002586-T1,OCOL_001021-T1,VNE69_01374-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 10.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002528-T1,COBT_001759-T1,GPK93_10g17610-T1,GINT2_001134-T1,GVAV_002610-T1,FUN_000444-T1,MDAP_001862-T1,OCOL_000206-T1,VNE69_07261-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 11.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001626-T1,COBT_000827-T1,GPK93_11g20630-T1,GINT2_000290-T1,GVAV_003434-T1,FUN_3529-T1,MDAP_000306-T1,OCOL_001203-T1,VNE69_11081-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 12.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001113-T1,COBT_001757-T1,GPK93_10g17830-T1,GINT2_001092-T1,GVAV_000067-T1,FUN_001200-T1,MDAP_002690-T1,OCOL_000184-T1,VNE69_07248-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 13.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001958-T1,COBT_000599-T1,GPK93_03g04040-T1,GINT2_002181-T1,GVAV_002467-T1,FUN_001753-T1,MDAP_000764-T1,OCOL_000330-T1,VNE69_10101-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 14.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000446-T1,COBT_000128-T1,GPK93_04g05970-T1,GINT2_000548-T1,GVAV_002219-T1,FUN_001401-T1,MDAP_002156-T1,OCOL_001036-T1,VNE69_04180-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 15.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002187-T1,COBT_000581-T1,GPK93_11g20040-T1,GINT2_000102-T1,GVAV_002493-T1,FUN_001979-T1,MDAP_001503-T1,OCOL_001158-T1,VNE69_01210-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 16.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000450-T1,COBT_002239-T1,GPK93_11g21110-T1,GINT2_001882-T1,GVAV_002207-T1,FUN_3534-T1,MDAP_002728-T1,OCOL_001245-T1,VNE69_01127-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 17.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001805-T1,COBT_002273-T1,GPK93_03g04520-T1,GINT2_001663-T1,GVAV_002022-T1,FUN_001293-T1,MDAP_002688-T1,OCOL_000242-T1,VNE69_01355-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 18.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001559-T1,COBT_003329-T1,GPK93_05g07280-T1,GINT2_001431-T1,GVAV_001010-T1,FUN_003044-T1,MDAP_000762-T1,OCOL_000689-T1,VNE69_09104-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 19.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001985-T1,COBT_000223-T1,GPK93_03g04290-T1,GINT2_000641-T1,GVAV_002935-T1,FUN_001819-T1,MDAP_001981-T1,OCOL_000351-T1,VNE69_10142-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 20.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002523-T1,COBT_000752-T1,GPK93_10g17720-T1,GINT2_000985-T1,GVAV_002601-T1,FUN_000429-T1,MDAP_002736-T1,OCOL_000195-T1,VNE69_07266-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 21.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001212-T1,COBT_000455-T1,GPK93_09g16390-T1,GINT2_001304-T1,GVAV_001137-T1,FUN_001505-T1,MDAP_001511-T1,OCOL_001478-T1,VNE69_08110-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 22.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001796-T1,COBT_002691-T1,GPK93_10g19380-T1,GINT2_002089-T1,GVAV_002034-T1,FUN_003254-T1,MDAP_000387-T1,OCOL_000043-T1,VNE69_08069-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 23.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001433-T1,COBT_003416-T1,GPK93_01g00210-T1,GINT2_001146-T1,GVAV_001082-T1,FUN_002882-T1,MDAP_002671-T1,OCOL_000488-T1,VNE69_06226-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 24.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000172-T1,COBT_000324-T1,GPK93_11g20650-T1,GINT2_000500-T1,GVAV_001312-T1,FUN_000683-T1,MDAP_000568-T1,OCOL_001205-T1,VNE69_01253-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 25.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000545-T1,COBT_001493-T1,GPK93_10g18600-T1,GINT2_001000-T1,GVAV_000631-T1,FUN_000882-T1,MDAP_000658-T1,OCOL_000117-T1,VNE69_07183-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 26.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000733-T1,COBT_000343-T1,GPK93_09g16470-T1,GINT2_000100-T1,GVAV_001896-T1,FUN_003284-T1,MDAP_000778-T1,OCOL_001486-T1,VNE69_08140-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 27.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002681-T1,COBT_000939-T1,GPK93_10g18390-T1,GINT2_000047-T1,GVAV_000577-T1,FUN_3422-T1,MDAP_002265-T1,OCOL_000134-T1,VNE69_07194-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 28.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000560-T1,COBT_001531-T1,GPK93_10g19010-T1,GINT2_000191-T1,GVAV_000657-T1,FUN_000248-T1,MDAP_000527-T1,OCOL_000074-T1,VNE69_07133-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 29.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000866-T1,COBT_000100-T1,GPK93_06g10400-T1,GINT2_001085-T1,GVAV_003369-T1,FUN_003207-T1,MDAP_000428-T1,OCOL_000749-T1,VNE69_05019-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 30.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001457-T1,COBT_001266-T1,GPK93_05g07750-T1,GINT2_000832-T1,GVAV_001114-T1,FUN_002077-T1,MDAP_002365-T1,OCOL_000645-T1,VNE69_05101-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 31.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002875-T1,COBT_003532-T1,GPK93_05g08010-T1,GINT2_001596-T1,GVAV_001225-T1,FUN_002224-T1,MDAP_002520-T1,OCOL_000622-T1,VNE69_05174-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 32.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002691-T1,COBT_003217-T1,GPK93_10g18570-T1,GINT2_001011-T1,GVAV_000588-T1,FUN_000284-T1,MDAP_002286-T1,OCOL_000119-T1,VNE69_07185-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 33.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001441-T1,COBT_000702-T1,GPK93_05g07730-T1,GINT2_000834-T1,GVAV_001087-T1,FUN_002066-T1,MDAP_001361-T1,OCOL_000647-T1,VNE69_05104-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 34.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001558-T1,COBT_000417-T1,GPK93_05g07270-T1,GINT2_001499-T1,GVAV_001009-T1,FUN_003042-T1,MDAP_000503-T1,OCOL_000690-T1,VNE69_09165-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 35.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001561-T1,COBT_003072-T1,GPK93_05g07290-T1,GINT2_000274-T1,GVAV_001026-T1,FUN_003049-T1,MDAP_002809-T1,OCOL_000688-T1,VNE69_09103-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 36.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002268-T1,COBT_001053-T1,GPK93_09g16000-T1,GINT2_001212-T1,GVAV_000881-T1,FUN_002962-T1,MDAP_000036-T1,OCOL_000390-T1,VNE69_01105-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 37.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001426-T1,COBT_003284-T1,GPK93_04g05700-T1,GINT2_002106-T1,GVAV_000305-T1,FUN_001155-T1,MDAP_002238-T1,OCOL_001013-T1,VNE69_01359-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 38.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001094-T1,COBT_000744-T1,GPK93_01g01160-T1,GINT2_001902-T1,GVAV_000125-T1,FUN_001235-T1,MDAP_000060-T1,OCOL_000571-T1,VNE69_02200-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 39.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001721-T1,COBT_000808-T1,GPK93_08g13240-T1,GINT2_001302-T1,GVAV_002670-T1,FUN_001515-T1,MDAP_000385-T1,OCOL_001665-T1,VNE69_02049-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 40.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001647-T1,COBT_002540-T1,GPK93_10g19180-T1,GINT2_000275-T1,GVAV_003473-T1,FUN_001438-T1,MDAP_001707-T1,OCOL_000061-T1,VNE69_09082-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 41.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002530-T1,COBT_000638-T1,GPK93_08g14940-T1,GINT2_001131-T1,GVAV_002612-T1,FUN_000448-T1,MDAP_002111-T1,OCOL_001804-T1,VNE69_07277-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 42.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001445-T1,COBT_003582-T1,GPK93_01g00150-T1,GINT2_000871-T1,GVAV_001098-T1,FUN_002867-T1,MDAP_001587-T1,OCOL_000482-T1,VNE69_06234-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 43.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001793-T1,COBT_003022-T1,GPK93_02g02420-T1,GINT2_001829-T1,GVAV_002032-T1,FUN_002600-T1,MDAP_001076-T1,OCOL_001317-T1,VNE69_12194-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 44.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002336-T1,COBT_001645-T1,GPK93_05g07860-T1,GINT2_000825-T1,GVAV_000770-T1,FUN_002921-T1,MDAP_000206-T1,OCOL_000635-T1,VNE69_05117-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 45.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002774-T1,COBT_003363-T1,GPK93_04g05510-T1,GINT2_000391-T1,GVAV_001470-T1,FUN_000534-T1,MDAP_000294-T1,OCOL_001001-T1,VNE69_01343-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 46.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001326-T1,COBT_001972-T1,GPK93_02g01900-T1,GINT2_000016-T1,GVAV_001045-T1,FUN_003080-T1,MDAP_001592-T1,OCOL_001361-T1,VNE69_06064-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 47.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001359-T1,COBT_000603-T1,GPK93_04g06340-T1,GINT2_001633-T1,GVAV_002143-T1,FUN_002515-T1,MDAP_001741-T1,OCOL_001069-T1,VNE69_12125-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 48.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001723-T1,COBT_000409-T1,GPK93_08g13150-T1,GINT2_001947-T1,GVAV_002667-T1,FUN_001507-T1,MDAP_002794-T1,OCOL_001657-T1,VNE69_02017-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 49.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001438-T1,COBT_001093-T1,GPK93_07g12210-T1,GINT2_001544-T1,GVAV_000251-T1,FUN_001333-T1,MDAP_001549-T1,OCOL_000869-T1,VNE69_04085-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 50.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000764-T1,COBT_001009-T1,GPK93_03g03740-T1,GINT2_000910-T1,GVAV_001873-T1,FUN_001663-T1,MDAP_000048-T1,OCOL_000304-T1,VNE69_11070-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 51.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001390-T1,COBT_001419-T1,GPK93_04g06090-T1,GINT2_001329-T1,GVAV_002111-T1,FUN_001775-T1,MDAP_001260-T1,OCOL_001046-T1,VNE69_01291-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 52.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001436-T1,COBT_002882-T1,GPK93_08g14460-T1,GINT2_000112-T1,GVAV_002926-T1,FUN_000283-T1,MDAP_002538-T1,OCOL_001766-T1,VNE69_07177-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 53.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001697-T1,COBT_003051-T1,GPK93_06g09990-T1,GINT2_000715-T1,GVAV_003158-T1,FUN_001593-T1,MDAP_002380-T1,OCOL_000783-T1,VNE69_02098-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 54.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000880-T1,COBT_000002-T1,GPK93_05g08340-T1,GINT2_001727-T1,GVAV_003357-T1,FUN_001162-T1,MDAP_000911-T1,OCOL_000039-T1,VNE69_05025-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 55.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001373-T1,COBT_000072-T1,GPK93_04g06260-T1,GINT2_001759-T1,GVAV_002125-T1,FUN_002527-T1,MDAP_002187-T1,OCOL_001062-T1,VNE69_12127-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 56.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002850-T1,COBT_000301-T1,GPK93_06g10010-T1,GINT2_000740-T1,GVAV_001173-T1,FUN_001611-T1,MDAP_000623-T1,OCOL_000781-T1,VNE69_02091-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 57.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002793-T1,COBT_001319-T1,GPK93_07g11230-T1,GINT2_000875-T1,GVAV_001440-T1,FUN_002649-T1,MDAP_002587-T1,OCOL_001626-T1,VNE69_03364-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 58.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000449-T1,COBT_003239-T1,GPK93_11g21100-T1,GINT2_000919-T1,GVAV_002208-T1,FUN_3536-T1,MDAP_000346-T1,OCOL_001244-T1,VNE69_01126-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 59.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001279-T1,COBT_003319-T1,GPK93_07g12490-T1,GINT2_002140-T1,GVAV_002204-T1,FUN_002583-T1,MDAP_001399-T1,OCOL_000895-T1,VNE69_01281-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 60.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001923-T1,COBT_000445-T1,GPK93_03g03810-T1,GINT2_002220-T1,GVAV_002868-T1,FUN_001697-T1,MDAP_001815-T1,OCOL_000311-T1,VNE69_11083-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 61.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000524-T1,COBT_001446-T1,GPK93_04g05720-T1,GINT2_000379-T1,GVAV_000581-T1,FUN_000501-T1,MDAP_000302-T1,OCOL_001015-T1,VNE69_01382-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 62.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000071-T1,COBT_001002-T1,GPK93_02g01800-T1,GINT2_002015-T1,GVAV_000412-T1,FUN_003055-T1,MDAP_002507-T1,OCOL_001370-T1,VNE69_06040-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 63.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000051-T1,COBT_003075-T1,GPK93_11g20460-T1,GINT2_001473-T1,GVAV_000424-T1,FUN_000797-T1,MDAP_002215-T1,OCOL_001192-T1,VNE69_01180-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 64.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002877-T1,COBT_002788-T1,GPK93_05g08030-T1,GINT2_001595-T1,GVAV_001226-T1,FUN_002223-T1,MDAP_002521-T1,OCOL_000621-T1,VNE69_05175-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 65.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001941-T1,COBT_003478-T1,GPK93_05g08720-T1,GINT2_001792-T1,GVAV_002998-T1,FUN_001054-T1,MDAP_001391-T1,OCOL_000011-T1,VNE69_12025-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 66.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000194-T1,COBT_000318-T1,GPK93_07g12440-T1,GINT2_001009-T1,GVAV_002003-T1,FUN_001377-T1,MDAP_002349-T1,OCOL_000891-T1,VNE69_04047-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 67.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000070-T1,COBT_001201-T1,GPK93_05g08060-T1,GINT2_001291-T1,GVAV_000245-T1,FUN_001588-T1,MDAP_002392-T1,OCOL_000618-T1,VNE69_05145-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 68.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001079-T1,COBT_001277-T1,GPK93_02g01850-T1,GINT2_000028-T1,GVAV_000188-T1,FUN_002934-T1,MDAP_002602-T1,OCOL_001366-T1,VNE69_06031-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 69.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002075-T1,COBT_000191-T1,GPK93_05g07350-T1,GINT2_001820-T1,GVAV_003522-T1,FUN_001841-T1,MDAP_000194-T1,OCOL_000683-T1,VNE69_05072-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 70.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000238-T1,COBT_001831-T1,GPK93_11g19710-T1,GINT2_000405-T1,GVAV_001371-T1,FUN_000632-T1,MDAP_001207-T1,OCOL_001129-T1,VNE69_01290-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 71.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001105-T1,COBT_002447-T1,GPK93_01g01250-T1,GINT2_001911-T1,GVAV_000074-T1,FUN_001209-T1,MDAP_000812-T1,OCOL_000579-T1,VNE69_02232-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 72.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002859-T1,COBT_003529-T1,GPK93_05g08140-T1,GINT2_000861-T1,GVAV_001182-T1,FUN_002093-T1,MDAP_001358-T1,OCOL_000611-T1,VNE69_05128-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 73.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001933-T1,COBT_000169-T1,GPK93_03g03970-T1,GINT2_000572-T1,GVAV_003006-T1,FUN_001750-T1,MDAP_002858-T1,OCOL_000325-T1,VNE69_10079-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 74.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001643-T1,COBT_000402-T1,GPK93_05g08450-T1,GINT2_001648-T1,GVAV_003468-T1,FUN_001927-T1,MDAP_000989-T1,OCOL_000032-T1,VNE69_08047-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 75.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002306-T1,COBT_001670-T1,GPK93_03g03550-T1,GINT2_000906-T1,GVAV_000833-T1,FUN_002268-T1,MDAP_000726-T1,OCOL_000288-T1,VNE69_11056-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 76.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001804-T1,COBT_000574-T1,GPK93_10g18670-T1,GINT2_000114-T1,GVAV_001997-T1,FUN_000313-T1,MDAP_000305-T1,OCOL_000110-T1,VNE69_07176-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 77.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000230-T1,COBT_002593-T1,GPK93_09g16110-T1,GINT2_001224-T1,GVAV_001364-T1,FUN_001586-T1,MDAP_002028-T1,OCOL_000379-T1,VNE69_03119-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 78.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001750-T1,COBT_002630-T1,GPK93_07g11300-T1,GINT2_001312-T1,GVAV_002088-T1,FUN_000948-T1,MDAP_002715-T1,OCOL_000798-T1,VNE69_02129-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 79.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002071-T1,COBT_001244-T1,GPK93_07g11480-T1,GINT2_001989-T1,GVAV_003518-T1,FUN_001997-T1,MDAP_000377-T1,OCOL_000815-T1,VNE69_03308-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 80.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001855-T1,COBT_000380-T1,GPK93_02g02890-T1,GINT2_001765-T1,GVAV_002095-T1,FUN_002405-T1,MDAP_000918-T1,OCOL_001544-T1,VNE69_12070-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 81.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002114-T1,COBT_002189-T1,GPK93_02g02070-T1,GINT2_001252-T1,GVAV_001300-T1,FUN_003260-T1,MDAP_000297-T1,OCOL_001347-T1,VNE69_06083-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 82.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002070-T1,COBT_002291-T1,GPK93_07g11460-T1,GINT2_000787-T1,GVAV_003516-T1,FUN_001991-T1,MDAP_000164-T1,OCOL_000813-T1,VNE69_03301-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 83.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001452-T1,COBT_000670-T1,GPK93_10g17700-T1,GINT2_001103-T1,GVAV_001478-T1,FUN_001501-T1,MDAP_001577-T1,OCOL_000197-T1,VNE69_07264-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 84.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001692-T1,COBT_000476-T1,GPK93_05g07680-T1,GINT2_000523-T1,GVAV_002721-T1,FUN_3587-T1,MDAP_000055-T1,OCOL_000652-T1,VNE69_05096-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 85.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001966-T1,COBT_003435-T1,GPK93_03g04460-T1,GINT2_000657-T1,GVAV_002480-T1,FUN_001785-T1,MDAP_001892-T1,OCOL_000367-T1,VNE69_01204-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 86.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000055-T1,COBT_001555-T1,GPK93_02g02030-T1,GINT2_001036-T1,GVAV_000420-T1,FUN_002828-T1,MDAP_001848-T1,OCOL_001351-T1,VNE69_06092-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 87.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001854-T1,COBT_000381-T1,GPK93_09g15570-T1,GINT2_000450-T1,GVAV_002094-T1,FUN_000875-T1,MDAP_001708-T1,OCOL_000955-T1,VNE69_03064-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 88.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000390-T1,COBT_000871-T1,GPK93_09g16990-T1,GINT2_001474-T1,GVAV_001656-T1,FUN_001024-T1,MDAP_001666-T1,OCOL_000432-T1,VNE69_03212-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 89.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001676-T1,COBT_002587-T1,GPK93_11g21380-T1,GINT2_002261-T1,GVAV_003567-T1,FUN_002749-T1,MDAP_000837-T1,OCOL_001268-T1,VNE69_02290-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 90.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000180-T1,COBT_003233-T1,GPK93_08g13270-T1,GINT2_001652-T1,GVAV_001125-T1,FUN_001171-T1,MDAP_001724-T1,OCOL_001668-T1,VNE69_08037-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 91.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001801-T1,COBT_000483-T1,GPK93_03g03830-T1,GINT2_001250-T1,GVAV_000195-T1,FUN_001756-T1,MDAP_002214-T1,OCOL_000313-T1,VNE69_11096-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 92.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002387-T1,COBT_002129-T1,GPK93_08g13660-T1,GINT2_000122-T1,GVAV_000697-T1,FUN_000144-T1,MDAP_001935-T1,OCOL_001697-T1,VNE69_07079-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 93.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001570-T1,COBT_002345-T1,GPK93_05g07110-T1,GINT2_001511-T1,GVAV_000988-T1,FUN_003010-T1,MDAP_000533-T1,OCOL_000702-T1,VNE69_09181-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 94.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001903-T1,COBT_002159-T1,GPK93_11g21230-T1,GINT2_000882-T1,GVAV_001783-T1,FUN_003324-T1,MDAP_000913-T1,OCOL_001256-T1,VNE69_01114-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 95.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001806-T1,COBT_001189-T1,GPK93_08g13090-T1,GINT2_001945-T1,GVAV_001327-T1,FUN_001473-T1,MDAP_001337-T1,OCOL_001652-T1,VNE69_02031-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 96.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000582-T1,COBT_000303-T1,GPK93_08g13520-T1,GINT2_000185-T1,GVAV_000680-T1,FUN_000175-T1,MDAP_001104-T1,OCOL_001684-T1,VNE69_07113-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 97.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001567-T1,COBT_001595-T1,GPK93_09g17030-T1,GINT2_001980-T1,GVAV_001005-T1,FUN_001032-T1,MDAP_000733-T1,OCOL_000436-T1,VNE69_03219-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 98.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002777-T1,COBT_003488-T1,GPK93_04g05480-T1,GINT2_001657-T1,GVAV_001465-T1,FUN_000541-T1,MDAP_001118-T1,OCOL_000998-T1,VNE69_01346-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 99.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000537-T1,COBT_001000-T1,GPK93_08g14120-T1,GINT2_002085-T1,GVAV_000731-T1,FUN_000090-T1,MDAP_001613-T1,OCOL_001737-T1,VNE69_01222-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 100.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001349-T1,COBT_002475-T1,GPK93_04g06420-T1,GINT2_001713-T1,GVAV_002133-T1,FUN_002511-T1,MDAP_001798-T1,OCOL_001076-T1,VNE69_11126-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 101.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002602-T1,COBT_000294-T1,GPK93_06g09090-T1,GINT2_000449-T1,GVAV_000523-T1,FUN_002595-T1,MDAP_000457-T1,OCOL_001447-T1,VNE69_04202-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 102.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000228-T1,COBT_002591-T1,GPK93_01g00340-T1,GINT2_001382-T1,GVAV_001362-T1,FUN_002902-T1,MDAP_001629-T1,OCOL_000500-T1,VNE69_06197-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 103.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000189-T1,COBT_000186-T1,GPK93_11g20840-T1,GINT2_000481-T1,GVAV_001325-T1,FUN_000440-T1,MDAP_002647-T1,OCOL_001219-T1,VNE69_01273-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 104.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001718-T1,COBT_000806-T1,GPK93_08g13110-T1,GINT2_001290-T1,GVAV_002672-T1,FUN_001476-T1,MDAP_001092-T1,OCOL_001654-T1,VNE69_02022-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 105.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002836-T1,COBT_001049-T1,GPK93_02g02810-T1,GINT2_000772-T1,GVAV_002837-T1,FUN_002451-T1,MDAP_001292-T1,OCOL_001552-T1,VNE69_12080-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 106.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000551-T1,COBT_000140-T1,GPK93_10g18820-T1,GINT2_000069-T1,GVAV_000634-T1,FUN_000321-T1,MDAP_000612-T1,OCOL_000096-T1,VNE69_07155-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 107.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000534-T1,COBT_002760-T1,GPK93_11g19620-T1,GINT2_001294-T1,GVAV_000602-T1,FUN_000625-T1,MDAP_001894-T1,OCOL_001123-T1,VNE69_01312-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 108.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001109-T1,COBT_000875-T1,GPK93_06g10450-T1,GINT2_000730-T1,GVAV_000066-T1,FUN_000323-T1,MDAP_001971-T1,OCOL_000745-T1,VNE69_05024-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 109.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002198-T1,COBT_003723-T1,GPK93_04g06180-T1,GINT2_001330-T1,GVAV_002115-T1,FUN_002543-T1,MDAP_001011-T1,OCOL_001054-T1,VNE69_04125-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 110.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000347-T1,COBT_000622-T1,GPK93_06g10520-T1,GINT2_001592-T1,GVAV_002311-T1,FUN_003216-T1,MDAP_000728-T1,OCOL_000739-T1,VNE69_05012-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 111.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002879-T1,COBT_002518-T1,GPK93_03g04180-T1,GINT2_000606-T1,GVAV_001228-T1,FUN_001835-T1,MDAP_002013-T1,OCOL_000340-T1,VNE69_10153-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 112.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002608-T1,COBT_002632-T1,GPK93_01g01410-T1,GINT2_000376-T1,GVAV_000492-T1,FUN_000253-T1,MDAP_000196-T1,OCOL_000595-T1,VNE69_02289-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 113.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000447-T1,COBT_000133-T1,GPK93_04g05860-T1,GINT2_002217-T1,GVAV_002216-T1,FUN_3537-T1,MDAP_002280-T1,OCOL_001025-T1,VNE69_04171-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 114.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000521-T1,COBT_000470-T1,GPK93_10g18740-T1,GINT2_000060-T1,GVAV_000591-T1,FUN_000302-T1,MDAP_000272-T1,OCOL_000103-T1,VNE69_07166-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 115.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000227-T1,COBT_001577-T1,GPK93_10g18620-T1,GINT2_000045-T1,GVAV_001357-T1,FUN_000402-T1,MDAP_002341-T1,OCOL_000115-T1,VNE69_07252-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 116.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002760-T1,COBT_002007-T1,GPK93_04g05390-T1,GINT2_000461-T1,GVAV_001482-T1,FUN_000511-T1,MDAP_000303-T1,OCOL_000989-T1,VNE69_01368-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 117.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002672-T1,COBT_000533-T1,GPK93_06g10030-T1,GINT2_000736-T1,GVAV_000569-T1,FUN_002265-T1,MDAP_002739-T1,OCOL_000779-T1,VNE69_02081-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 118.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000049-T1,COBT_001208-T1,GPK93_02g03010-T1,GINT2_001781-T1,GVAV_000435-T1,FUN_000376-T1,MDAP_000210-T1,OCOL_001532-T1,VNE69_12059-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 119.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001028-T1,COBT_003440-T1,GPK93_09g16680-T1,GINT2_001388-T1,GVAV_000227-T1,FUN_003293-T1,MDAP_002873-T1,OCOL_001504-T1,VNE69_08164-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 120.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001219-T1,COBT_001007-T1,GPK93_03g03460-T1,GINT2_000929-T1,GVAV_003137-T1,FUN_003003-T1,MDAP_001783-T1,OCOL_000280-T1,VNE69_11033-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 121.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001859-T1,COBT_001708-T1,GPK93_07g12380-T1,GINT2_001534-T1,GVAV_000704-T1,FUN_001154-T1,MDAP_000257-T1,OCOL_000885-T1,VNE69_04023-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 122.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002078-T1,COBT_000195-T1,GPK93_07g11490-T1,GINT2_001991-T1,GVAV_003525-T1,FUN_002001-T1,MDAP_002554-T1,OCOL_000816-T1,VNE69_03307-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 123.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000482-T1,COBT_001895-T1,GPK93_08g14200-T1,GINT2_000734-T1,GVAV_001503-T1,FUN_000876-T1,MDAP_000672-T1,OCOL_001744-T1,VNE69_02018-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 124.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002762-T1,COBT_002628-T1,GPK93_04g05370-T1,GINT2_000459-T1,GVAV_001479-T1,FUN_000518-T1,MDAP_000352-T1,OCOL_000987-T1,VNE69_01365-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 125.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002779-T1,COBT_001920-T1,GPK93_05g08250-T1,GINT2_000050-T1,GVAV_001468-T1,FUN_001503-T1,MDAP_002190-T1,OCOL_000602-T1,VNE69_05134-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 126.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000210-T1,COBT_001983-T1,GPK93_10g18270-T1,GINT2_000077-T1,GVAV_000779-T1,FUN_000329-T1,MDAP_000799-T1,OCOL_000145-T1,VNE69_07186-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 127.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002223-T1,COBT_000262-T1,GPK93_06g09270-T1,GINT2_000847-T1,GVAV_003277-T1,FUN_001049-T1,MDAP_001305-T1,OCOL_001432-T1,VNE69_05141-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 128.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002192-T1,COBT_002729-T1,GPK93_05g08240-T1,GINT2_000866-T1,GVAV_001186-T1,FUN_001420-T1,MDAP_000639-T1,OCOL_000603-T1,VNE69_07033-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 129.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002274-T1,COBT_002913-T1,GPK93_01g00610-T1,GINT2_001230-T1,GVAV_001575-T1,FUN_002972-T1,MDAP_002872-T1,OCOL_000525-T1,VNE69_06160-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 130.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001367-T1,COBT_000030-T1,GPK93_04g05850-T1,GINT2_001445-T1,GVAV_002139-T1,FUN_001481-T1,MDAP_002498-T1,OCOL_001024-T1,VNE69_04170-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 131.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001018-T1,COBT_003780-T1,GPK93_02g02290-T1,GINT2_001174-T1,GVAV_000243-T1,FUN_003182-T1,MDAP_002351-T1,OCOL_001327-T1,VNE69_08022-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 132.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002761-T1,COBT_003361-T1,GPK93_04g05240-T1,GINT2_000463-T1,GVAV_001480-T1,FUN_000510-T1,MDAP_000168-T1,OCOL_000975-T1,VNE69_01369-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 133.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001017-T1,COBT_003560-T1,GPK93_07g12120-T1,GINT2_001557-T1,GVAV_000244-T1,FUN_3571-T1,MDAP_000549-T1,OCOL_000861-T1,VNE69_04096-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 134.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002880-T1,COBT_002941-T1,GPK93_05g08050-T1,GINT2_001603-T1,GVAV_001232-T1,FUN_002233-T1,MDAP_002547-T1,OCOL_000619-T1,VNE69_05180-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 135.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000184-T1,COBT_000189-T1,GPK93_03g03720-T1,GINT2_000911-T1,GVAV_001323-T1,FUN_001811-T1,MDAP_001846-T1,OCOL_000302-T1,VNE69_11076-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 136.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002178-T1,COBT_000636-T1,GPK93_09g17300-T1,GINT2_001370-T1,GVAV_003070-T1,FUN_001099-T1,MDAP_001772-T1,OCOL_000457-T1,VNE69_03240-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 137.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000448-T1,COBT_003327-T1,GPK93_11g21090-T1,GINT2_002207-T1,GVAV_002211-T1,FUN_3535-T1,MDAP_002090-T1,OCOL_001243-T1,VNE69_01129-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 138.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002594-T1,COBT_000041-T1,GPK93_08g13350-T1,GINT2_001369-T1,GVAV_000051-T1,FUN_002897-T1,MDAP_002598-T1,OCOL_001674-T1,VNE69_03242-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 139.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000033-T1,COBT_001634-T1,GPK93_01g01150-T1,GINT2_001904-T1,GVAV_000439-T1,FUN_002658-T1,MDAP_001059-T1,OCOL_000570-T1,VNE69_02199-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 140.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001637-T1,COBT_003032-T1,GPK93_04g06410-T1,GINT2_000802-T1,GVAV_001301-T1,FUN_001164-T1,MDAP_002785-T1,OCOL_001075-T1,VNE69_11127-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 141.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002062-T1,COBT_004141-T1,GPK93_04g06220-T1,GINT2_001956-T1,GVAV_003503-T1,FUN_001960-T1,MDAP_002724-T1,OCOL_001058-T1,VNE69_12135-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 142.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001710-T1,COBT_000282-T1,GPK93_05g07380-T1,GINT2_001806-T1,GVAV_001387-T1,FUN_003058-T1,MDAP_001793-T1,OCOL_000680-T1,VNE69_05039-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 143.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002773-T1,COBT_000675-T1,GPK93_04g05310-T1,GINT2_000455-T1,GVAV_001472-T1,FUN_002185-T1,MDAP_001093-T1,OCOL_000982-T1,VNE69_01337-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 144.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002534-T1,COBT_002172-T1,GPK93_08g14830-T1,GINT2_001114-T1,GVAV_002622-T1,FUN_000475-T1,MDAP_001628-T1,OCOL_001796-T1,VNE69_07292-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 145.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002847-T1,COBT_000302-T1,GPK93_07g11800-T1,GINT2_000869-T1,GVAV_001174-T1,FUN_000454-T1,MDAP_002516-T1,OCOL_000839-T1,VNE69_10007-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 146.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002506-T1,COBT_001833-T1,GPK93_10g17920-T1,GINT2_001081-T1,GVAV_002583-T1,FUN_002674-T1,MDAP_002325-T1,OCOL_000177-T1,VNE69_07255-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 147.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001986-T1,COBT_001683-T1,GPK93_03g04390-T1,GINT2_000636-T1,GVAV_002931-T1,FUN_001814-T1,MDAP_002376-T1,OCOL_000360-T1,VNE69_10129-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 148.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000730-T1,COBT_003444-T1,GPK93_02g02590-T1,GINT2_000402-T1,GVAV_000782-T1,FUN_003128-T1,MDAP_000074-T1,OCOL_001569-T1,VNE69_12154-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 149.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002219-T1,COBT_000802-T1,GPK93_02g02480-T1,GINT2_001823-T1,GVAV_001370-T1,FUN_002264-T1,MDAP_000731-T1,OCOL_001311-T1,VNE69_03113-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 150.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002610-T1,COBT_002048-T1,GPK93_08g12950-T1,GINT2_001147-T1,GVAV_003560-T1,FUN_001428-T1,MDAP_001319-T1,OCOL_001638-T1,VNE69_02005-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 151.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_000870-T1,COBT_003659-T1,GPK93_06g10460-T1,GINT2_001560-T1,GVAV_003365-T1,FUN_003189-T1,MDAP_000425-T1,OCOL_000744-T1,VNE69_06112-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 152.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002718-T1,COBT_002255-T1,GPK93_05g07530-T1,GINT2_000537-T1,GVAV_001056-T1,FUN_003096-T1,MDAP_000695-T1,OCOL_000665-T1,VNE69_05075-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 153.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_002838-T1,COBT_001412-T1,GPK93_11g20720-T1,GINT2_000238-T1,GVAV_002851-T1,FUN_002964-T1,MDAP_002584-T1,OCOL_001210-T1,VNE69_01245-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 154.fasta
proteinortho_grab_proteins.pl -exact 'BDAP_001943-T1,COBT_000171-T1,GPK93_03g03990-T1,GINT2_000880-T1,GVAV_002990-T1,FUN_003063-T1,MDAP_002825-T1,OCOL_000326-T1,VNE69_10109-T1' Binucleata_daphniae.faa Conglomerata_obtusa.faa Encephalitozoon_intestinalis.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_tvaerminnensis.faa Mitosporidium_daphniae.faa Ordospora_colligata.faa Vairimorpha_necatrix.faa > 155.fasta

#sorted.species.list
#>Binucleata daphniae
#>Conglomerata obtusa
#>Hamiltosporidium tvaerminnensis
#>Glugoides intestinalis
#>Encephalitozoon intestinalis
#>Gurleya vavrai
#>Mitosporidium daphniae
#>Ordospora colligata
#>Vairimorpha necatrix

# sort sequences by header and replace header with species name
for fasta in *.fasta
do
paste -d'\n' sorted.species.list <(awk 'NR%2==0' <(bioawk -c fastx '{print}' $fasta | sort -k1,1V | awk '{print ">"$1;print $2}')) > $fasta.sorted.fasta
done

# infer gene trees from SCOs (align, trim, make tree)
for fasta in *sorted.fasta
do
mafft --anysymbol --quiet $fasta > $fasta.mafft.fasta
sed -i "" 's/ /_/g' $fasta.mafft.fasta
trimal -in $fasta.mafft.fasta -out $fasta.mafft.trimal.phylip -automated1 -phylip
trimal -in $fasta.mafft.fasta -out $fasta.mafft.trimal.fasta -automated1 -fasta
iqtree2 -s $fasta.mafft.trimal.phylip -nt AUTO -ntmax 12 -seed 12345 -bb 1000 -o "Mitosporidium_daphniae"
done

# ASTRAL "weigthed" phylogeny based on gene trees
cat *fasta.sorted.fasta.mafft.trimal.phylip.treefile > all.fasta.sorted.fasta.mafft.trimal.phylip.treefile
~/bioinformatics/ASTER-MacOS/bin/wastral -r 16 -s 16 -x 100 -n 0 --root "Rozella_allomycis" -t 12 -o all.fasta.sorted.fasta.mafft.trimal.phylip.astral all.fasta.sorted.fasta.mafft.trimal.phylip.treefile 2>all.fasta.sorted.fasta.mafft.trimal.phylip.log

# concatenate all aligned and trimmed SCOs
#seqkit concat *sorted.fasta | sed -e 's/\|.*$//' > sorted.concat.fasta
cut -f1 *trimal.fasta | awk -v RS=">" -v FS="\n" -v OFS="\n" '{for(i=2; i<=NF; i++) {seq[$1] = seq[$1]$i}}; END {for(id in seq){print ">"id, seq[id]}}' > sorted.mafft.trimal.concat.fasta

# infer phylogeny
iqtree2 -s sorted.mafft.trimal.concat.fasta -nt AUTO -ntmax 12 -seed 12345 -bb 1000 -o "Mitosporidium_daphniae"
