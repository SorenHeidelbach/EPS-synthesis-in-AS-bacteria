#Script for PSI-BLAST search with operon as query

data="/srv/PHN/Projects/bt8cf21/data"
WD="/srv/PHN/Projects/bt8cf21/SorenH"
database=/srv/PHN/Projects/bt8cf21/MortenE/database/databasefolder
cd $WD
date=$(date '+%d_%b')
mkdir $WD/"psiblast_results"
mkdir $WD/"psiblast_results"/$date

##Query path definitions, non-MSA queries
operon_fasta=(
Operon_Alginate_algD-algA.fasta  
Operon_cellulose.fasta                  
Operon_Psl_PslA-PslL.fasta
Operon_Pel_PelA-PelG.fasta      
Operon_PNAG_Achromobacter_aegrifaciens.fasta  
Operon_PNAG_x5.fasta    
Operon_PNAG_E_coli_NP_415540-43.1.fasta       
Operon_PNAG_Xanthomonas_vasicola.fasta  
Operon_PNAG_Variovorax_sp_RO1.fasta           
Operon_PNAG_Y_Pestis_AAB66588-91.fasta  
Operon_succinoglycan.fasta
Operon_Xanthan.fasta)

## PSI-BLAST of all operons in the GFF database from HQ-MAGs
module load BLAST+/2.11.0-foss-2020b
for operon in ${operon_fasta[@]}; do
psiblast -query $WD/"query_input/single_genes/"$operon -db $database -out $WD/psiblast_results/$date/$operon -evalue 0.0001 -qcov_hsp_perc 25 -max_hsps 1 -max_target_seqs 5000 -outfmt 6 -num_iterations 3 -comp_based_stats 1 -num_threads 5
done
module purge




