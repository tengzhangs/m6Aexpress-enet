conda activate R4
cd ./sramp_simple/
##MOCKL
nohup perl runsramp.pl ./SRAMP/full_RNA_mode/exomePeak2_result/new_motif_peak_seq.fa ./exomePeak2_result/singlebase_m6Asites.txt full &
