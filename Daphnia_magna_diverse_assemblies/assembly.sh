## example

## assembly with hifiasm (v.0.19.8-r603)

hifiasm -o DZ-OER-1.l0.asm -l0 -t8 ../Pool_DE_Circular_Consensus_Sequencing_Reads/m84177_240307_192125_s3.bc2082--bc2082.hifi_reads.fastq.gz

## cleaning with ncbi-fcs (v.0.5.0)

# database info:
#gx-version: build:Jan 19 2023 15:50:20; git:v0.3.0-151-g9aad15db
#build-date: Tue Jan 24 16:18:22 EST 2023

source ~/software/fcs_env/bin/activate # python v.3.8
export FCS_DEFAULT_IMAGE=~/software/ncbi-fcs/fcs-gx.sif 
python3 ~/software/ncbi-fcs/fcs.py screen genome --fasta DZ-OER-1.asm.bp.p_ctg.fa --out-dir ./gx_out/ --gx-db ~/software/ncbi-fcs/gxdb --tax-id 35525
cat DZ-OER-1.asm.bp.p_ctg.fa | python3 ~/software/ncbi-fcs/fcs.py clean genome --action-report ./gx_out/DZ-OER-1.asm.bp.p_ctg.35525.fcs_gx_report.txt --output DZ-OER-1.asm.bp.p_ctg.clean.fasta --contam-fasta-out DZ-OER-1.asm.bp.p_ctg.contam.fasta


