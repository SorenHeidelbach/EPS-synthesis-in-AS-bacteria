

DB="/srv/PHN/Projects/bt8cf21/data/MGP1000_HQMAG1083_prot_db_prokv1.14_20200210.faa"
SPLIT_DB="/srv/PHN/Projects/bt8cf21/MGP1000_HQMAG1083_prot_db_split"
mkdir $SPLIT_DB
cat $DB | awk -v SPLIT_DB=$SPLIT_DB '/^>/ {MAG=sprintf(SPLIT_DB"/%s.fasta",substr($0,2));
split(MAG, a, " ");
OUT = substr(a[1], 1, length(a[1])-6)};
OUT {print > OUT".faa"}'

# When interpreting the script think of going through each line.
# First it is checked if the line start with >, if it does, then the line is saved to the variable MAG and processed to remove prokka name (a[1]) and id (OUT)
# The WHOLE line is then printed to the file "OUT" (variable with only MAG id)
# Then the next line is evaluated for ">" and printed to the file "OUT". As the aa sequence that follow the ">" indentifiers, is the aa sequence, it is printed to the same file below the ">" identifier
# [The "OUT {print > OUT".faa"}" means print the whole line always]
# This is done for each line until a new > is found, which now overrides the old MAG, a[1] and OUT variables

