cut -f1 pulexOpsins_vs_CLR4Genome.tblastn  | sort | uniq -c | sed 's/^ *//g'| sed 's/ /,/g' | awk -F, ' ($1 < 3000) ' | cut -f2 -d "," > good_tblastn_id ; for i in `cat good_tblastn_id` ; do grep "^$i	" pulexOpsins_vs_CLR4Genome.tblastn >> Filtered_pulexOpsins_vs_CLR4Genome.tblastn ; done 
mv Filtered_pulexOpsins_vs_CLR4Genome.tblastn pulexOpsins_vs_CLR4Genome.tblastn

