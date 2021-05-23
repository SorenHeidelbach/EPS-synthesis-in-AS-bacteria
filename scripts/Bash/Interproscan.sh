
# Working directionary, our group folder
WD=/srv/PHN/Projects/bt8cf21
cd $WD
# Folder with all subsetted fasta files
FASTA=/srv/PHN/Projects/bt8cf21/fasta_output
# Results folder
RESULTS=$WD/results/interproscan_Pfam

# Loading interproscan
module load InterProScan/5.38-76.0-foss-2018a

for j in $FASTA/*; do
# For each polysaccharide
j=`basename $j`
echo $j
mkdir -p $RESULTS/$j
mkdir -p $RESULTS/$j\_fasta_removed
for i in $FASTA/$j/*; do
	# For each MAG id
	# Removing .faa from the name
	i=`basename $i | rev | cut -c 5- | rev`
	echo $i
	if [[ ! -s $RESULTS/$j\_fasta_removed/$i.gff3 ]]; then
	# Runnign interproscan
	## Only Pfam and SUPERFAMILY database
	## Only GFF3 file is written out
	interproscan.sh -cpu 24 -appl Pfam -f GFF3 -o $RESULTS/$j/$i.gff3 -i $FASTA/$j/$i.faa
	
	# Remove trailing fasta sequence in gff3 file and leading three lines
	F=$RESULTS/$j/$i.gff3
	sed '/^##FASTA$/,$d' ${F} | sed '1,/^##interproscan*/d' > $RESULTS/$j\_fasta_removed/$i.gff3
	fi
done
done
