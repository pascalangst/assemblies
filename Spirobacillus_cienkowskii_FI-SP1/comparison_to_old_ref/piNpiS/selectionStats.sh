# Identify orthologs
proteinortho6.pl Spirobacillus_cienkowskii_FI-SP1.faa comparison_to_old_ref/GCA_003339605.1/protein.faa -project=comparison_to_old_ref/proteinortho

# cut orthologs from Genbank assembly
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < GCA_003339605.1/cds_from_genomic.fna > cds_old.oneline.fna
sed -i "" -E 's/(>)lcl.Q............._cds_(..........)[^ ]* .*/\1\2|KSL3/g' cds_old.oneline.fna 
for i in `awk '{print $2}' proteinortho/proteinortho.single-copy.tsv`;do awk 'BEGIN{RS=">"} /'$i'/{print ">" $0}' cds_old.oneline.fna;done > cds_old.orthologs.fasta
paste <(awk '{print $2}' proteinortho/proteinortho.single-copy.tsv) <(awk '{print $1}' proteinortho/proteinortho.single-copy.tsv) | while read a b; do sed -i "" "s/$a/$b/" cds_old.orthologs.fasta; done

# cut orthologs from new assembly
sed -E 's/(>Spiro2_.....).*/\1|FI-SP1/g' ../Spirobacillus_cienkowskii_FI-SP1.ffn > cds_new.oneline.fna
for i in `awk '{print $1}' proteinortho/proteinortho.single-copy.tsv`;do awk 'BEGIN{RS=">"} /'$i'/{print ">" $0}' cds_new.oneline.fna;done > cds_new.orthologs.fasta

# create one fasta for each assembly and ortholog
awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' cds_old.orthologs.fasta
mkdir fasta_references
for f in *Spiro2*fasta; do mv "$f" "$(echo "fasta_references/$f" | sed -e 's/>//;s/.fasta//')";done

awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' cds_old.orthologs.fasta
mkdir fasta_references/fasta_old
for f in *Spiro2*fasta; do mv "$f" "$(echo "fasta_references/fasta_old/$f" | sed -e 's/>//')";done

awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' cds_new.orthologs.fasta
mkdir fasta_references/fasta_new
for f in *Spiro2*fasta; do mv "$f" "$(echo "fasta_references/fasta_new/$f" | sed -e 's/>//')";done

# align sequences (these are the input for selectionStats.py)
ls Spiro2_* | parallel -j 6 "Rscript ../align.R {}"

# concatenate all aligned fasta
for sample in *aligned.fas; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $sample > $sample.oneline.fasta; done
paste *aligned.fas.oneline.fasta > all.fasta
sed -E "s/    //g" all.fasta > all2.fasta
sed -E 's/(>[^>]*)>.*/\1/g' all2.fasta > all.fasta
sed -i 1d all.fasta
sed -i "" 1d all.fasta

# calculate hamming distance
python ../Hammingdistance.py 
#0.011845542556730525