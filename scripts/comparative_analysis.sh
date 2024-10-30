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



cd /home/pascal/assemblies/comp_geno_PhD

conda activate mamba-funannotate

export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/pasa-2.4.1

funannotate compare --cpus 12 --input ~/Gi_seq_second/funannotate_IL-YERU-16/annotate_results/Glugoides_intestinalis_IL-YERU-16.gbk \
~/assemblies/Mito_SP/masurca_ont_pb/CA.mr.99.17.15.0.02/contigs.SP.braker.dikarya/ \
~/assemblies/genebank/Paramicrosporidium_saccamoebae/ \
~/assemblies/genebank/Spraguea_lophii/ \
~/assemblies/genebank/Anncaliia_algerae/ \
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
~/assemblies/genebank/Hamiltosporidium_magnivora/ \
~/assemblies/genebank/Tubulinosema_ratisbonensis/ \
~/assemblies/genebank/Dictyocoela_muelleri/ \
~/assemblies/genebank/Thelohania_contejeani/ \
~/assemblies/genebank/Cucumispora_dikerogammari/ \
~/BE-OM-3_Bd-2-2_depl/CA.mr.49.17.15.0.02/funannotate_BE-OM-3.v1/ \
~/Gv_SP/nextDenovo_polished_sgs-lgs_medaka/funannotate_Gv_SP.v2/ \
~/assemblies/SK-39_Larssonia/megahit/funannotate_FI-SK-39.v1.1/ \
~/assemblies/hifi/GB-LK1-1/manual_merge/funannotate_GB-LK1-1/ \
~/Ht_methylation_project/analysis/iso-seq/long_read_protocol/FI-OER-3-3.v2.6/
#-------------------------------------------------------
#[Aug 13 11:16 AM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.12
#[Aug 13 11:16 AM]: Running 1.8.11
#[Aug 13 11:16 AM]: Now parsing 28 genomes
#[Aug 13 11:16 AM]: working on Glugoides intestinalis
#[Aug 13 11:16 AM]: working on Mitosporidium daphniae
#[Aug 13 11:16 AM]: working on Paramicrosporidium saccamoebae
#[Aug 13 11:16 AM]: working on Spraguea lophii 42_110
#[Aug 13 11:16 AM]: working on Anncaliia algerae PRA339
#[Aug 13 11:16 AM]: working on Encephalitozoon romaleae SJ-2008
#[Aug 13 11:16 AM]: working on Nematocida displodere
#[Aug 13 11:16 AM]: working on Nematocida parisii ERTm1
#[Aug 13 11:16 AM]: working on Nosema ceranae
#[Aug 13 11:16 AM]: working on Rozella allomycis CSF55
#[Aug 13 11:16 AM]: working on Vavraia culicis subsp. floridensis
#[Aug 13 11:17 AM]: working on Vittaforma corneae ATCC 50505
#[Aug 13 11:17 AM]: working on Edhazardia aedis USNM 41457
#[Aug 13 11:17 AM]: working on Nosema bombycis CQ1
#[Aug 13 11:17 AM]: working on Enterocytozoon bieneusi H348
#[Aug 13 11:17 AM]: working on Hepatospora eriocheir
#[Aug 13 11:17 AM]: working on Trachipleistophora hominis
#[Aug 13 11:17 AM]: working on Pseudoloma neurophilia
#[Aug 13 11:17 AM]: working on Hamiltosporidium magnivora
#[Aug 13 11:17 AM]: working on Tubulinosema ratisbonensis
#[Aug 13 11:17 AM]: working on Dictyocoela muelleri
#[Aug 13 11:17 AM]: working on Thelohania contejeani
#[Aug 13 11:17 AM]: working on Cucumispora dikerogammari
#[Aug 13 11:18 AM]: working on Binucleata daphniae
#[Aug 13 11:18 AM]: working on Gurleya daphniae
#[Aug 13 11:18 AM]: working on Conglomerata obtusa
#[Aug 13 11:18 AM]: working on Ordospora colligata
#[Aug 13 11:18 AM]: working on Hamiltosporidium tvaerminnensis
#[Aug 13 11:18 AM]: No secondary metabolite annotations found
#[Aug 13 11:18 AM]: Summarizing PFAM domain results
#Fontconfig warning: ignoring UTF-8: not a valid region tag
#[Aug 13 11:18 AM]: Summarizing InterProScan results
#[Aug 13 11:18 AM]: Loading InterPro descriptions
#[Aug 13 11:18 AM]: Summarizing MEROPS protease results
#[Aug 13 11:18 AM]: found 16/87 MEROPS familes with stdev >= 1.000000
#[Aug 13 11:18 AM]: Summarizing CAZyme results
#[Aug 13 11:18 AM]: found 4/56 CAZy familes with stdev >= 1.000000
#[Aug 13 11:18 AM]: Summarizing COG results
#[Aug 13 11:18 AM]: Summarizing secreted protein results
#[Aug 13 11:18 AM]: Summarizing fungal transcription factors
#[Aug 13 11:18 AM]: Running GO enrichment for each genome
#  WARNING: skipping Hamiltosporidium_magnivora.txt as no GO terms
#  WARNING: skipping Edhazardia_aedis_USNM_41457.txt as no GO terms
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
#[Aug 13 11:19 AM]: Running orthologous clustering tool, ProteinOrtho.  This may take awhile...
#[Aug 13 11:28 AM]: Compiling all annotations for each genome
#[Aug 13 11:29 AM]: Inferring phylogeny using iqtree
#[Aug 13 11:29 AM]: Found 6 single copy BUSCO orthologs, will use all to infer phylogeny
#[Aug 13 12:41 PM]: Compressing results to output file: funannotate_compare.tar.gz
#[Aug 13 12:43 PM]: Funannotate compare completed successfully!

# do phylogeny based on the identified SCOs
proteinortho_grab_proteins.pl -exact 'H312_00699-T1,BDAP_001326-T1,LOBT_001972-T1,CDIK_1973-T1,DMUE_2428-T1,EDEG_03183-T1,EROM_020250-T1,EBI_21688-T1,GINT_000359-T1,GVAV_001045-T1,CWI36_0678p0010-T1,FUN_003080-T1,HERIO_1499-T1,MDAP_001592-T1,NEDG_00370-T1,NEPG_00197-T1,NBO_4g0069-T1,G9O61_00g018800-T1,OCOL_001361-T1,PSACC_00668-T1,M153_13560001237-T1,ROZALSC1DRAFT_29250-T1,SLOPH_7-T1,TCON_0101-T1,THOM_2344-T1,TUBRATIS_22680-T1,VCUG_01595-T1,VICG_00714-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 1.fasta
proteinortho_grab_proteins.pl -exact 'H312_00789-T1,BDAP_002877-T1,LOBT_002788-T1,CDIK_3368-T1,DMUE_0402-T1,EDEG_01293-T1,EROM_050950-T1,EBI_27165-T1,GINT_001875-T1,GVAV_001226-T1,CWI36_0280p0020-T1,FUN_002223-T1,HERIO_1310-T1,MDAP_002521-T1,NEDG_00128-T1,NEPG_00519-T1,NBO_72g0004-T1,G9O61_00g016100-T1,OCOL_000621-T1,PSACC_01191-T1,M153_3560004869-T1,ROZALSC1DRAFT_27522-T1,SLOPH_48-T1,TCON_2510-T1,THOM_2953-T1,TUBRATIS_13470-T1,VCUG_02140-T1,VICG_01910-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 2.fasta
proteinortho_grab_proteins.pl -exact 'H312_02923-T1,BDAP_000450-T1,LOBT_002239-T1,CDIK_1788-T1,DMUE_5550-T1,EDEG_02940-T1,EROM_111410-T1,EBI_24670-T1,GINT_002253-T1,GVAV_002207-T1,CWI36_0599p0050-T1,FUN_3534-T1,HERIO_924-T1,MDAP_002728-T1,NEDG_00805-T1,NEPG_01108-T1,NBO_10g0045-T1,G9O61_00g002220-T1,OCOL_001245-T1,PSACC_01566-T1,M153_6500024225-T1,ROZALSC1DRAFT_28270-T1,SLOPH_1999-T1,TCON_2475-T1,THOM_3138-T1,TUBRATIS_001240-T1,VCUG_01638-T1,VICG_01056-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 3.fasta
proteinortho_grab_proteins.pl -exact 'H312_03446-T1,BDAP_001639-T1,LOBT_002683-T1,CDIK_1959-T1,DMUE_4985-T1,EDEG_03668-T1,EROM_060660-T1,EBI_23921-T1,GINT_000028-T1,GVAV_003455-T1,CWI36_1572p0010-T1,FUN_002694-T1,HERIO_2697-T1,MDAP_001529-T1,NEDG_02068-T1,NEPG_02208-T1,NBO_6g0063-T1,G9O61_00g002920-T1,OCOL_001408-T1,PSACC_02884-T1,M153_19030001338-T1,ROZALSC1DRAFT_28835-T1,SLOPH_1562-T1,TCON_0514-T1,THOM_2610-T1,TUBRATIS_001680-T1,VCUG_00451-T1,VICG_00759-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 4.fasta
proteinortho_grab_proteins.pl -exact 'H312_01630-T1,BDAP_002774-T1,LOBT_003363-T1,CDIK_1150-T1,DMUE_2113-T1,EDEG_01226-T1,EROM_040400-T1,EBI_22833-T1,GINT_000678-T1,GVAV_001470-T1,CWI36_3044p0010-T1,FUN_000534-T1,HERIO_1446-T1,MDAP_000294-T1,NEDG_00316-T1,NEPG_00291-T1,NBO_554g0006-T1,G9O61_00g011750-T1,OCOL_001001-T1,PSACC_03723-T1,M153_1100042559-T1,ROZALSC1DRAFT_27770-T1,SLOPH_25-T1,TCON_1062-T1,THOM_1010-T1,TUBRATIS_000190-T1,VCUG_00901-T1,VICG_00544-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 5.fasta
proteinortho_grab_proteins.pl -exact 'H312_03235-T1,BDAP_001452-T1,LOBT_000670-T1,CDIK_0066-T1,DMUE_2713-T1,EDEG_01612-T1,EROM_100090-T1,EBI_27253-T1,GINT_000903-T1,GVAV_001478-T1,CWI36_0560p0040-T1,FUN_001501-T1,HERIO_1929-T1,MDAP_001577-T1,NEDG_01648-T1,NEPG_02004-T1,NBO_507g0006-T1,G9O61_00g004630-T1,OCOL_000197-T1,PSACC_02324-T1,M153_1364000159-T1,ROZALSC1DRAFT_29961-T1,SLOPH_2456-T1,TCON_1274-T1,THOM_1444-T1,TUBRATIS_16190-T1,VCUG_00041-T1,VICG_01546-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 6.fasta
proteinortho_grab_proteins.pl -exact 'H312_01305-T1,BDAP_001359-T1,LOBT_000603-T1,CDIK_2702-T1,DMUE_1538-T1,EDEG_01201-T1,EROM_041160-T1,EBI_27095-T1,GINT_001834-T1,GVAV_002143-T1,CWI36_1861p0020-T1,FUN_002515-T1,HERIO_1181-T1,MDAP_001741-T1,NEDG_00621-T1,NEPG_00865-T1,NBO_32g0045-T1,G9O61_00g009880-T1,OCOL_001069-T1,PSACC_00724-T1,M153_1570003290-T1,ROZALSC1DRAFT_28387-T1,SLOPH_345-T1,TCON_0250-T1,THOM_1775-T1,TUBRATIS_003730-T1,VCUG_00672-T1,VICG_00792-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 7.fasta
proteinortho_grab_proteins.pl -exact 'H312_02492-T1,BDAP_001373-T1,LOBT_000072-T1,CDIK_0632-T1,DMUE_5184-T1,EDEG_03076-T1,EROM_041080-T1,EBI_24216-T1,GINT_001716-T1,GVAV_002125-T1,CWI36_0268p0050-T1,FUN_002527-T1,HERIO_1178-T1,MDAP_002187-T1,NEDG_00614-T1,NEPG_00857-T1,NBO_32g0048-T1,G9O61_00g009840-T1,OCOL_001062-T1,PSACC_01998-T1,M153_7190003529-T1,ROZALSC1DRAFT_30462-T1,SLOPH_1979-T1,TCON_0249-T1,THOM_0419-T1,TUBRATIS_003770-T1,VCUG_02110-T1,VICG_00784-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 8.fasta
proteinortho_grab_proteins.pl -exact 'H312_00419-T1,BDAP_000210-T1,LOBT_001983-T1,CDIK_0550-T1,DMUE_1251-T1,EDEG_02786-T1,EROM_100630-T1,EBI_22821-T1,GINT_000294-T1,GVAV_000779-T1,CWI36_0568p0030-T1,FUN_000329-T1,HERIO_2138-T1,MDAP_000799-T1,NEDG_01364-T1,NEPG_01782-T1,NBO_375g0006-T1,G9O61_00g005330-T1,OCOL_000145-T1,PSACC_00101-T1,M153_13690002889-T1,ROZALSC1DRAFT_27300-T1,SLOPH_2461-T1,TCON_0158-T1,THOM_2146-T1,TUBRATIS_002770-T1,VCUG_00546-T1,VICG_00073-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 9.fasta
proteinortho_grab_proteins.pl -exact 'H312_01347-T1,BDAP_000733-T1,LOBT_000343-T1,CDIK_0117-T1,DMUE_3920-T1,EDEG_00570-T1,EROM_090840-T1,EBI_23131-T1,GINT_000269-T1,GVAV_001896-T1,CWI36_1968p0010-T1,FUN_003284-T1,HERIO_2427-T1,MDAP_000778-T1,NEDG_00705-T1,NEPG_00976-T1,NBO_41g0021-T1,G9O61_00g000870-T1,OCOL_001486-T1,PSACC_00193-T1,M153_2750007406-T1,ROZALSC1DRAFT_18919-T1,SLOPH_1446-T1,TCON_0839-T1,THOM_2306-T1,TUBRATIS_23360-T1,VCUG_01454-T1,VICG_00203-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 10.fasta
proteinortho_grab_proteins.pl -exact 'H312_02256-T1,BDAP_002336-T1,LOBT_001645-T1,CDIK_2229-T1,DMUE_4804-T1,EDEG_00716-T1,EROM_050790-T1,EBI_24873-T1,GINT_001086-T1,GVAV_000770-T1,CWI36_0215p0030-T1,FUN_002921-T1,HERIO_521-T1,MDAP_000206-T1,NEDG_01458-T1,NEPG_01244-T1,NBO_53g0029-T1,G9O61_00g009230-T1,OCOL_000635-T1,PSACC_02963-T1,M153_2200020223-T1,ROZALSC1DRAFT_31056-T1,SLOPH_310-T1,TCON_0317-T1,THOM_1380-T1,TUBRATIS_008340-T1,VCUG_00209-T1,VICG_01603-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 11.fasta
proteinortho_grab_proteins.pl -exact 'H312_03238-T1,BDAP_001643-T1,LOBT_000402-T1,CDIK_0806-T1,DMUE_0650-T1,EDEG_01378-T1,EROM_051340-T1,EBI_23231-T1,GINT_000708-T1,GVAV_003468-T1,CWI36_0114p0040-T1,FUN_001927-T1,HERIO_27-T1,MDAP_000989-T1,NEDG_00898-T1,NEPG_02477-T1,NBO_609gi001-T1,G9O61_00g000360-T1,OCOL_000032-T1,PSACC_00251-T1,M153_4127000461-T1,ROZALSC1DRAFT_29142-T1,SLOPH_2294-T1,TCON_2685-T1,THOM_3095-T1,TUBRATIS_23980-T1,VCUG_01731-T1,VICG_00597-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 12.fasta
proteinortho_grab_proteins.pl -exact 'H312_00444-T1,BDAP_002682-T1,LOBT_000940-T1,CDIK_1276-T1,DMUE_2053-T1,EDEG_03573-T1,EROM_100710-T1,EBI_26123-T1,GINT_000319-T1,GVAV_000578-T1,CWI36_0328p0040-T1,FUN_000275-T1,HERIO_1880-T1,MDAP_000158-T1,NEDG_01329-T1,NEPG_01739-T1,NBO_84g0001-T1,G9O61_00g005280-T1,OCOL_000138-T1,PSACC_00497-T1,M153_1890003809-T1,ROZALSC1DRAFT_29469-T1,SLOPH_2478-T1,TCON_0145-T1,THOM_2264-T1,TUBRATIS_12430-T1,VCUG_01152-T1,VICG_00081-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 13.fasta
proteinortho_grab_proteins.pl -exact 'H312_01166-T1,BDAP_000524-T1,LOBT_001446-T1,CDIK_0267-T1,DMUE_2627-T1,EDEG_01657-T1,EROM_040550-T1,EBI_22378-T1,GINT_000690-T1,GVAV_000581-T1,CWI36_1222p0020-T1,FUN_000501-T1,HERIO_61-T1,MDAP_000302-T1,NEDG_01152-T1,NEPG_01329-T1,NBO_64g0025-T1,G9O61_00g012100-T1,OCOL_001015-T1,PSACC_03732-T1,M153_2500035249-T1,ROZALSC1DRAFT_27593-T1,SLOPH_1701-T1,TCON_2572-T1,THOM_0227-T1,TUBRATIS_003500-T1,VCUG_00940-T1,VICG_00509-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 14.fasta
proteinortho_grab_proteins.pl -exact 'H312_02108-T1,BDAP_001437-T1,LOBT_000551-T1,CDIK_0540-T1,DMUE_2398-T1,EDEG_02744-T1,EROM_081510-T1,EBI_21958-T1,GINT_002388-T1,GVAV_000137-T1,CWI36_0156p0060-T1,FUN_000113-T1,HERIO_789-T1,MDAP_002125-T1,NEDG_01888-T1,NEPG_01597-T1,NBO_12g0011-T1,G9O61_00g022470-T1,OCOL_001771-T1,PSACC_02817-T1,M153_7060003767-T1,ROZALSC1DRAFT_30734-T1,SLOPH_2047-T1,TCON_0242-T1,THOM_1367-T1,TUBRATIS_16960-T1,VCUG_00195-T1,VICG_00889-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 15.fasta
proteinortho_grab_proteins.pl -exact 'H312_00416-T1,BDAP_001806-T1,LOBT_001189-T1,CDIK_0516-T1,DMUE_4879-T1,EDEG_00441-T1,EROM_080160-T1,EBI_24859-T1,GINT_002204-T1,GVAV_001327-T1,CWI36_2446p0020-T1,FUN_001473-T1,HERIO_382-T1,MDAP_001337-T1,NEDG_01757-T1,NEPG_01907-T1,NBO_59g0003-T1,G9O61_00g017520-T1,OCOL_001652-T1,PSACC_03357-T1,M153_27110001801-T1,ROZALSC1DRAFT_6479-T1,SLOPH_2200-T1,TCON_0245-T1,THOM_1918-T1,TUBRATIS_16600-T1,VCUG_02635-T1,VICG_01804-T1' Anncaliia_algerae_PRA339.faa Binucleata_daphniae.faa Conglomerata_obtusa.faa Cucumispora_dikerogammari.faa Dictyocoela_muelleri.faa Edhazardia_aedis_USNM_41457.faa Encephalitozoon_romaleae_SJ-2008.faa Enterocytozoon_bieneusi_H348.faa Glugoides_intestinalis.faa Gurleya_daphniae.faa Hamiltosporidium_magnivora.faa Hamiltosporidium_tvaerminnensis.faa Hepatospora_eriocheir.faa Mitosporidium_daphniae.faa Nematocida_displodere.faa Nematocida_parisii_ERTm1.faa Nosema_bombycis_CQ1.faa Nosema_ceranae.faa Ordospora_colligata.faa Paramicrosporidium_saccamoebae.faa Pseudoloma_neurophilia.faa Rozella_allomycis_CSF55.faa Spraguea_lophii_42_110.faa Thelohania_contejeani.faa Trachipleistophora_hominis.faa Tubulinosema_ratisbonensis.faa Vavraia_culicis_subsp._floridensis.faa Vittaforma_corneae_ATCC_50505.faa > 16.fasta




#sorted.species.list
#>Cucumispora dikerogammari
#>Hamiltosporidium magnivora
#>Dictyocoela muelleri
#>Enterocytozoon bieneusi
#>Edhazardia aedis
#>Encephalitozoon romaleae
#>Hamiltosporidium tvaerminnensis
#>Nosema ceranae
#>Glugoides intestinalis
#>Gurleya vavrai
#>Anncaliia algerae
#>Hepatospora_eriocheir
#>Conglomerata obtusa
#>Pseudoloma_neurophilia
#>Mitosporidium daphniae
#>Nosema bombycis
#>Nematocida displodere
#>Nematocida_parisii
#>Ordospora colligata
#>Paramicrosporidium saccamoebae
#>Rozella allomycis
#>Spraguea lophii
#>Astathelohania contejeani
#>Trachipleistophora hominis
#>Tubulinosema ratisbonensis
#>Vavraia culicis subsp. floridensis
#>Vittaforma corneae


# sort sequences by header and replace header with species name
for fasta in *.fasta
do
paste -d'\n' sorted.species.list <(awk 'NR%2==0' <(bioawk -c fastx '{print}' $fasta | sort -k1,1V | awk '{print ">"$1;print $2}')) > $fasta.sorted.fasta
done

# concatenate all SCOs and align them
#seqkit concat *sorted.fasta | sed -e 's/\|.*$//' > sorted.concat.fasta
cut -f1 -d "-" *.sorted.fasta | awk -v RS=">" -v FS="\n" -v OFS="\n" '{for(i=2; i<=NF; i++) {seq[$1] = seq[$1]$i}}; END {for(id in seq){print ">"id, seq[id]}}' > sorted.concat.fasta
mafft --anysymbol --quiet sorted.concat.fasta > sorted.concat.mafft.fasta

# trim alignment
sed -i 's/ /_/g' sorted.concat.mafft.fasta
trimal -in sorted.concat.mafft.fasta -out sorted.concat.trimal.phylip -automated1 -phylip

# infer phylogeny
iqtree2 -s sorted.concat.trimal.phylip -nt AUTO -ntmax 12 -seed 12345 -bb 1000 -o "Rozella_allomycis"
