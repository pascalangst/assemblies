curl -LO https://github.com/ncbi/fcs/raw/main/dist/fcs.py
curl https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/latest/fcs-gx.sif -Lo fcs-gx.sif
export FCS_DEFAULT_IMAGE=~/software/ncbi-fcs/fcs-gx.sif 
python3 fcs.py db get --mft https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest --dir gxdb
cat AUS.asm.bp.p_ctg.fa | python3 ~/software/ncbi-fcs/fcs.py clean genome --action-report ./gx_out/AUS.asm.bp.p_ctg.85468.fcs_gx_report.txt --output AUS.asm.bp.p_ctg.clean.fasta --contam-fasta-out AUS.asm.bp.p_ctg.contam.fasta
