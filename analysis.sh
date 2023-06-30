
#!/bin/bash

# Analysis script, currently does not work

# Access the passed variables
outdir="$1"
fastqs="$2"
slanes="$3"
prot="$4"
pval="$5"
qval="$6"
tval="$7"
threads="$8"
extra="$9"
lval="${10}"
dval="${11}"

# move to working directory
# make output directories
cd $outdir
mkdir counts individual.lanes fastqs fastas histos 2>/dev/null

# loop through reads & process them
for R1 in $fastqs/*R1*
do

	# Define some general use variables
	basename=$(basename ${R1})
	lbase=${basename//_R*}
	sbase=${basename//_L00*}
	R2=${R1//R1_001.fastq*/R2_001.fastq*}

	# Make a 'sample' directory for all analyses
	# and combined lane outputs (aka 'sample' outputs)
	dir=$outdir/individual.lanes/$sbase
	mkdir $dir 2>/dev/null

	 # Make a directory for indiv lane read & histo outputs
	lhist=$dir/histos
	fadir=$dir/fastas
	fqdir=$dir/fastqs
	cdir=$dir/counts
	mkdir $lhist $fadir $fqdir $cdir 2>/dev/null

	# Join reads & extract insert
	echo "Joining $lbase reads & extracting primer..."
	pandaseq -f $R1 -r $R2 -F \
	$pval $qval \
	-w $fqdir/$lbase.joined.fastq $tval -T $threads $extra $lval $dval 2>/dev/null

	# Convert to fasta
	echo "Converting joined $lbase FASTQ to FASTA..."
	awk 'NR%4 == 1 {print ">" substr($0, 2)} NR%4 == 2 {print}' $fqdir/$lbase.joined.fastq > $fadir/$lbase.joined.fasta

	# Combine sequences from each lane into single files
	echo "Adding $lbase reads to total $sbase reads..."
	cat $fqdir/$lbase.joined.fastq >> $dir/$sbase.joined.fastq
    cat $fadir/$lbase.joined.fasta >> $dir/$sbase.joined.fasta

	# Generate indiv lane  nt length distributions
	if [ -z $slanes ];
        then
        	:
		else
			echo "Generating $lbase nt length distribution for individual lanes..."
			# Length distribution for nt sequences
			awk 'NR%2 == 0 {reads+=1}END{print "#Reads:", reads}' $fadir/$lbase.joined.fasta > $lhist/$lbase.joined.nt.histo
			awk 'NR == 2 {line = $0; min = length($1)} NR > 2 && NR%2 ==0 && length($1) < min {line = $0; min = length($1)} END{print "Min:", min}' $fadir/$lbase.joined.fasta >> $lhist/$lbase.joined.nt.histo
			awk 'NR == 2 {line = $0; counts = length($1)} NR > 2 && NR%2 ==0 && length($1) > max {line = $0; max = length($1)} END{print "Max:", max}' $fadir/$lbase.joined.fasta >> $lhist/$lbase.joined.nt.histo

			# Read Length Histogram:
			echo "#Read Length Histogram:" >> $lhist/$lbase.joined.nt.histo
			echo "Len" "Reads" "%Reads" | column -t >> $lhist/$lbase.joined.nt.histo
			awk 'NR%2 == 0  {reads+=1; a[length($1)]+=1}END{for (i in a) printf "%s %s %.3f%% \n",i, a[i], 100*a[i]/reads}' $fadir/$lbase.joined.fasta | column -t | sort -g -k1  >> $lhist/$lbase.joined.nt.histo

			# Counts for nt lanes
			echo "Calculating unique & total reads for lane $lbase..."
			unique=$(awk 'NR%2 == 0 {seen[$1] += 1} END {print length(seen)}' $fadir/$lbase.joined.fasta)
			total=$(awk 'NR%2 == 0 {tot += 1} END {print tot}' $fadir/$lbase.joined.fasta)

			# Collect unique, total and sequences with counts in the counts file
			echo "Collecting unique, total and sequences in file..."
			echo "number of unique sequences = $unique" > $cdir/$lbase.joined.counts.txt
			echo "total number of molecules = $total"   >> $cdir/$lbase.joined.counts.txt
			echo  >>  $cdir/$lbase.joined.counts.txt
			awk -v tot="$total" 'NR%2 == 0 {seen[$1] += 1} END {for (i in seen) printf "%s %s %.3f%% \n", i, seen[i], 100*seen[i]/tot}' $fadir/$lbase.joined.fasta| column -t | sort -n -r -k2 >>  $cdir/$lbase.joined.counts.txt
			echo ""
	fi
done

cd $outdir/individual.lanes

########## CREATE COUNTS FILE FOR DNA ##########

# Loop through directories and generate count files

ls -1 | while read d
do
	test -d "$d" || continue
	# Define variables
	base=$(basename ${d})
	(cd $d ;

	# Calculate unique and total
	echo "Calculating unique & total reads for $base..."
	unique=$(awk 'NR%2 == 0 {seen[$1] += 1} END {print length(seen)}' $base.joined.fasta)
	total=$(awk 'NR%2 == 0 {tot += 1} END {print tot}' $base.joined.fasta)

	# Collect unique, total and sequences with counts in the counts file
	echo "number of unique sequences = $unique" > $outdir/counts/$base\_counts.txt
	echo "total number of molecules = $total"   >> $outdir/counts/$base\_counts.txt
	echo  >>  $outdir/counts/$base\_counts.txt
	awk -v tot="$total" 'NR%2 == 0 {seen[$1] += 1} END {for (i in seen) printf "%s %s %.3f%% \n", i, seen[i], 100*seen[i]/tot}' $base.joined.fasta | column -t | sort -n -r -k2 >>  $outdir/counts/$base\_counts.txt

	# Redirect outputs
	mv $base.joined.fasta $outdir/fastas/$base.joined.fasta
	mv $base.joined.fastq $outdir/fastqs/$base.joined.fastq

    )
done

# Cleanup indiv lanes

echo ""

if [ -z $slanes ];
	then
		echo "Cleaning up all individual lane outputs..."
                rm -r $outdir/individual.lanes/
        else
               	echo "Individual lane outputs will be retained"
fi

########## CREATE HISTO FILE FOR DNA ##########

cd $outdir/counts

for file in *counts.txt

do
	# Generate DNA length distribution for all lanes combined
	echo ""
	echo "Generating ${file//_counts.txt} DNA length distribution..."
	awk 'NR > 3 {reads+=$2}END{print "#Reads:", reads}' $file > ${file//_counts.txt}'_counts_histo.txt'
	awk 'NR == 4 {line = $0; min = length($1)} NR > 4 && length($1) < min {line = $0; min = length($1)} END{print "Min:", min}' $file >> ${file//_counts.txt}'_counts_histo.txt'
	awk 'NR == 4 {line = $0; max = length($1)} NR > 4 && length($1) > max {line = $0; max = length($1)} END{print "Max:", max}' $file >> ${file//_counts.txt}'_counts_histo.txt'

	# Read Length Histogram:
	echo "#Read Length Histogram:" >> ${file//_counts.txt}'_counts_histo.txt'
	echo "Len" "Reads" "%Reads" | column -t >> ${file//_counts.txt}'_counts_histo.txt'
	awk 'NR > 3  {reads+=$2; a[length($1)]+=$2}END{for (i in a) printf "%s %s %.3f%% \n",i, a[i], 100*a[i]/reads}' $file | column -t | sort -g -k1  >> ${file//_counts.txt}'_counts_histo.txt'

	mv ${file//_counts.txt}'_counts_histo.txt' ../histos/${file//_counts.txt}'_counts_histo.txt'

########## TRANSLATE DNA INTO PEPTIDES ##########

if [ -z $prot ];
	then
		echo "No translation will be performed"

	else

	# Translate into aa
	echo "Translating ${file//_counts.txt} DNA to peptides..."
	# temp=$(pwd) # Save our current working directory
	# cd "$hdir"
	python  ../../../translator.py $file # "$(echo "$hdir/translator.py")" "$(echo "$temp/$file")" # TODO: Change this to be less scuffed? It breaks when data is in a different file
	# cd "$temp"

	# Print in new file every line except the first 3 (2 with the number of molecules and sequences and town empty lines):
	tail -n +4 ${file//_counts.txt}'_counts.aa.dup.txt' | sort > newfile.txt;

	#Calculate total reads
	awk '{totaa+=$2}END{print totaa}' ${file//_counts.txt}'_counts.aa.dup.txt'  > /dev/null

	# Remove duplicates and sum abundances:
	awk '{seen[$1]+=$2;totaa+=$2}END{for (i in seen) printf  "%s %s %.3f%%\n", i, seen[i], 100*seen[i]/totaa}' newfile.txt  |column -t | sort -k1 > newfile2.txt;

	# Sort following abundance:
	sort newfile2.txt -k2 -n -r > newfile3.txt;

	# Print second line:
	val=$(awk 'BEGIN {FS = " "} ; {sum+=$2} END {print sum}' newfile3.txt) ; echo total number of molecules      =     $val | cat - newfile3.txt > newfile4.txt;

	# Print first line:
	val=$(cat newfile3.txt | wc -l) ; echo number of unique sequences     =     $val | cat - newfile4.txt > newfile5.txt;

	# Print third (empty) line:
	awk 'NR==3 {print ""} 1' newfile5.txt > ${file//_counts.txt}'_counts.aa.txt';

	# Remove every temp file:
	rm newfile.txt; rm newfile2.txt; rm newfile3.txt; rm newfile4.txt; rm newfile5.txt
	rm ${file//_counts.txt}'_counts.aa.dup.txt'

########## CREATE HISTO FILE FOR PEPTIDES ##########

	# Generate peptide length distribution for all lanes combined
	echo "Generating ${file//_counts.txt} aa length distribution..."
	awk 'NR > 3 {reads+=$2}END{print "#Reads:", reads}' ${file//_counts.txt}'_counts.aa.txt' > ${file//_counts.txt}'_counts.aa_histo.txt'
	awk 'NR == 4 {line = $0; min = length($1)} NR > 4 && length($1) < min {line = $0; min = length($1)} END{print "Min:", min}' ${file//_counts.txt}'_counts.aa.txt' >> ${file//_counts.txt}'_counts.aa_histo.txt'
	awk 'NR == 4 {line = $0; max = length($1)} NR > 4 && length($1) > max {line = $0; max = length($1)} END{print "Max:", max}' ${file//_counts.txt}'_counts.aa.txt' >> ${file//_counts.txt}'_counts.aa_histo.txt'

	# Read Length Histogram:
	echo "#Read Length Histogram:" >>  ${file//_counts.txt}'_counts.aa_histo.txt'
	echo "Len" "Reads" "%Reads" | column -t >>  ${file//_counts.txt}'_counts.aa_histo.txt'
	awk 'NR > 3  {reads+=$2; a[length($1)]+=$2}END{for (i in a) printf "%s %s %.3f%% \n",i, a[i], 100*a[i]/reads}' ${file//_counts.txt}'_counts.aa.txt' | column -t | sort -g -k1  >>  ${file//_counts.txt}'_counts.aa_histo.txt'

fi

done

cd ..

if [ -z $prot ];
	then
		:
	else
		mkdir counts.aa
		mv counts/*aa.txt counts.aa/
		mv counts/*aa.dup.txt counts.aa/
		mv counts/*aa_histo.txt histos/
fi

########## CREATE LOG FILE FOR DNA ##########

if [ -z $prot ];

	then

		cd ..

		echo ""  >> $outdir/log.txt
		echo "sample" "fastq_R1" "fastq_R2" "unique_nt" "total_nt" "recovered_nt(%)"| column -t > $outdir/log_temp1.txt

		for R1 in *R1*
		do
			basename=$(basename ${R1})
			lbase=${basename//_R*}
			sbase=${basename//_L00*}
			R2=${R1//R1_001.fastq*/R2_001.fastq*}

			echo $sbase \
			$(cat $R1 | zcat | awk 'END {print NR/4}') \
			$(cat $R2 | zcat | awk 'END {print NR/4}')  \
			$(cat ${outdir}/counts/$sbase\_counts.txt | awk 'BEGIN {ORS=" "}; NR==1{print $6}' ) \
			$(cat ${outdir}/counts/$sbase\_counts.txt | awk 'BEGIN {ORS=" "}; NR==2{print $6}' ) \
			| column -t >> $outdir/log_temp2.txt
		done

		awk '{ printf "%s %.2f%%\n", $0, 100*$5/$2 }' $outdir/log_temp2.txt | column -t >> $outdir/log_temp1.txt
		awk  '{print }' $outdir/log_temp1.txt | column -t >> $outdir/log.txt

		rm $outdir/log_temp1.txt
		rm $outdir/log_temp2.txt

	else


########## CREATE LOG FILE FOR DNA AND PEPTIDE SEQUENCES ##########

		cd ..

		echo ""  >> $outdir/log.txt
		echo "sample" "fastq_R1" "fastq_R2" "unique_nt" "total_nt" "recovered_nt(%)" "unique_aa" "total_aa" "recovered_aa(%)"| column -t > $outdir/log_temp1.txt

		for R1 in *R1*
		do

			basename=$(basename ${R1})
			lbase=${basename//_R*}
			sbase=${basename//_L00*}
			R2=${R1//R1_001.fastq*/R2_001.fastq*}

			echo $sbase \
			$(cat $R1 | zcat | awk 'END {print NR/4}') \
			$(cat $R2 | zcat | awk 'END {print NR/4}')  \
			$(cat ${outdir}/counts/$sbase\_counts.txt | awk 'BEGIN {ORS=" "}; NR==1{print $6}' ) \
			$(cat ${outdir}/counts/$sbase\_counts.txt | awk 'BEGIN {ORS=" "}; NR==2{print $6}' ) \
			$(cat ${outdir}/counts.aa/$sbase\_counts.aa.txt | awk 'BEGIN {ORS=" "}; NR==1{print $6}' ) \
			$(cat ${outdir}/counts.aa/$sbase\_counts.aa.txt | awk 'BEGIN {ORS=" "}; NR==2{print $6}' ) \
			| column -t >> $outdir/log_temp2.txt

		done

		awk '{ printf "%s %s %s %s %s %.2f%% %s %s %.2f%%\n", $1, $2, $3, $4, $5, 100*$5/$2, $6, $7, 100*$7/$2 }' $outdir/log_temp2.txt | column -t >> $outdir/log_temp1.txt
		awk  '{print }' $outdir/log_temp1.txt | column -t >> $outdir/log.txt

		rm $outdir/log_temp1.txt
		rm $outdir/log_temp2.txt

fi

