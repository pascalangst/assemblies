sed -i "" 's/(REVERSE SENSE)/_reverse/g' orf_list.fasta
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
