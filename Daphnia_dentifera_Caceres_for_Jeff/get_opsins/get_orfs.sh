# get pulex opsins
seqtk subseq ncbi_dataset/ncbi_dataset/data/GCA_000187875.1/protein.faa pulex_opsins.list > pulex_opsins.fa

# tblastn with pulex opsin genes against the genome, with an evalue of 1e-5
makeblastdb -in CLR4.asm.bp.p_ctg.clean.fasta -parse_seqids -dbtype nucl -out CLR4_db
tblastn -query pulex_opsins.fa -db CLR4_db -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out pulexOpsins_vs_CLR4Genome.tblastn -num_threads 6
cut -f1 pulexOpsins_vs_CLR4Genome.tblastn  | sort | uniq -c | sed 's/^ *//g'| sed 's/ /,/g' | awk -F, ' ($1 < 3000) ' | cut -f2 -d "," > good_tblastn_id ; for i in `cat good_tblastn_id` ; do grep "^$i	" pulexOpsins_vs_CLR4Genome.tblastn >> Filtered_pulexOpsins_vs_CLR4Genome.tblastn ; done 
mv Filtered_pulexOpsins_vs_CLR4Genome.tblastn pulexOpsins_vs_CLR4Genome.tblastn

# Extract non-overlapping best-hit sequences (first round of tblastn with the evalue of 1e-5) and extend 1000bp upstream and 1000bp downstream
Rscript extract_bestHit_and_extend.R pulexOpsins_vs_CLR4Genome.tblastn 
samtools faidx CLR4.asm.bp.p_ctg.clean.fasta < Best_hits_filtered.tsv > Best_hits_filtered.fasta
sed -i "" 's/:/-/g' Best_hits_filtered.fasta

# Extract ORFs
getorf -sequence Best_hits_filtered.fasta -outseq orf_list.fasta -minsize 750 -find 3

#Rename results of getorf in order to have good fasta headers

sed -i 's/(REVERSE SENSE)/_reverse/g' orf_list.fasta
grep ">" orf_list.fasta | sed 's/>//g' > oldfastaheaders
sed 's/-/	/g' oldfastaheaders | sed 's/_[0-9] \[/	/g' | sed 's/\]//g' | sed 's/_reverse/	reverse/g'  > ren_oldfastaheaders

IFS=$'\n' #treat the file line by line in the following for loop

for line in `cat ren_oldfastaheaders` ; do
		scaffold=`echo "$line" | cut -f1`

		if grep -q "reverse" <<< "$line" ; then 
			coord_start=`echo "$line" | awk 'BEGIN {FS="	"} {sum+=$2+$5-1} END {print sum}'`
			coord_end=`echo "$line" | awk 'BEGIN {FS="	"} {sum+=$2+$4-1} END {print sum}'`
		else
			coord_start=`echo "$line" | awk 'BEGIN {FS="	"} {sum+=$2+$4-1} END {print sum}'`
			coord_end=`echo "$line" | awk 'BEGIN {FS="	"} {sum+=$2+$5-1} END {print sum}'`

		fi

		echo "$scaffold-$coord_start-$coord_end" 


done > newfastaheaders

paste -d "\t" oldfastaheaders newfastaheaders > renaming_file

perl rename_fasta.pl renaming_file orf_list.fasta > renamed_orf_list.fasta

#remove identical sequences
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' renamed_orf_list.fasta | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > renamed_orf_list.fasta_uniq.fa