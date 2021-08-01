# Uses .tsv file genrated from "psiblast_subset" variable in main R script and the split aa database

WD="/srv/PHN/Projects/bt8cf21/SorenH/"
cd $WD
PSI="./psiblast_subset_reduced.tsv"
DB="/srv/PHN/Projects/bt8cf21/MGP1000_HQMAG1083_prot_db_split"


# Generate fasta files by unique operon
DIR="/srv/PHN/Projects/bt8cf21/SorenH/AA_operon"
mkdir $DIR
OPERONS=`cat $PSI | awk '{print $5}' | sort | uniq -d`
for i in ${OPERONS[@]}; do
MAGS=`cat $PSI | awk -v i=$i '{if ($5 == i){print $3}}' | sort | uniq`
for j in ${MAGS[@]}; do
	PROKKA=`cat $PSI | awk -v j=$j -v i=$i '{if ($3 == j && $5 == i) {print $2}}' | sort | uniq`
	cat $DB/$j.faa | while read line; do
	
	if [[ $line == *">"* ]]; then
	IFS=" "; read -a magid <<< "$line"; # Split line by " "
	magid="${magid[0]}"; # Remove proposed gene
	magid="${magid:1}" # remove ">"
	GENE=`cat $PSI | awk -v magid=$magid '{if ($2 == magid) {print $1}}' | sort | uniq`
	MAG_TAX=`cat $PSI | awk -v magid=$magid '{if ($2 == magid) {print "p."$7"_f."$10"_g."$11}}' | uniq`
	line=">$GENE $MAG_TAX $magid"
	fi
	
	if [[ " ${PROKKA[@]} " =~ "$magid" ]]; then
	echo $line >> $DIR/$j.operon_$i.fa
	fi
	done
done
done


# Generate fasta files by gene name
DIR="/srv/PHN/Projects/bt8cf21/SorenH/AA_gene"
mkdir $DIR
GENES=`cat $PSI | awk '{print $1}' | sort | uniq -d`
for i in ${GENES[@]}; do
MAG=`cat $PSI | awk -v i=$i '{if ($1 == i){print $3}}' | sort | uniq`
for j in ${MAG[@]}; do
	PROKKA=`cat $PSI | awk -v j=$j -v i=$i '{if ($3 == j && $1 == i) {print $2}}' | sort | uniq`
	cat $DB/$j.faa | while read line; do
	
	# Identify sequence header -> extract MAG+PROKKA to magid -> generate new sequence header with MAG+PROKKA+GENE
	if [[ $line == *">"* ]]; then
	IFS=" "; read -a magid <<< "$line"; magid="${magid[0]}";	magid="${magid:1}"
	GENE=`cat $PSI | awk -v magid=$magid '{if ($2 == magid) {print $1}}' | sort | uniq`
	line=">$magid $GENE"
	fi
	
	# If "magid" in array of MAG+PROKKA, print line to fasta file
	if [[ " ${PROKKA[@]} " =~ "$magid" ]]; then
	echo $line >> $DIR/gene_$i.fa
	fi
	done
done
done


# Generate fasta files by unique operon and the surrounding genes
DIR="/srv/PHN/Projects/bt8cf21/SorenH/AA_operon_with_surrounding_genes"
mkdir $DIR
OPERONS=`cat $PSI | awk '{print $5}' | sort | uniq -d`
for i in ${OPERONS[@]}; do
MAGS=`cat $PSI | awk -v i=$i '{if ($5 == i){print $3}}' | sort | uniq`
for j in ${MAGS[@]}; do
	MAG_PROKKA=`cat $PSI | awk -v j=$j -v i=$i '{if ($3 == j && $5 == i) {print $2}}' | sort | uniq`
	PROKKA=`cat $PSI | awk -v j=$j -v i=$i '{if ($3 == j && $5 == i) {print $4}}' | sort | uniq`
	# Get max Prokka ID for operon
	IFS=$'\n'
	max=`echo "${PROKKA[*]}" | sort -nr | head -n1`
	let "max = $max + 1"
	
	# Get min Prokka ID for operon
	min=`echo "${PROKKA[*]}" | sort -nr | tail -n1`
	let "min = $min + 1"
	
	cat $DB/$j.faa | while read line; do
	if [[ $line == *">"* ]]; then
	IFS=" "; read -a magid <<< "$line" # Split line by " "
	PROKKA_annotation="${magid[@]:1}"
	magid="${magid[0]}" # Remove proposed gene
	magid="${magid:1}" # remove ">"
	PROKKA_id_line="${magid:(-5)}" 
	PROKKA_id_line=`echo $PROKKA_id_line | sed "s/^0*\([1-9]\)/\1/;s/^0*$/0/"`
	# extract only prokka id
	
	GENE=`cat $PSI | awk -v magid=$magid '{if ($2 == magid) {print $1}}' | sort | uniq`
	MAG_TAX=`cat $PSI | awk -v magid=$magid '{if ($2 == magid) {print "p."$7"_f."$10"_g."$11}}' | uniq`
	line=">$GENE $MAG_TAX $magid $PROKKA_annotation"
	fi
	
	if [[ $PROKKA_id_line -le $max && $PROKKA_id_line -ge $min ]]; then
	echo $line >> $DIR/$j.operon_$i.fa
	fi
	done
done
done


interproscan.sh -appl Pfam -b ./interproscan/ -i ./AA_operon_with_surrounding_genes/*
