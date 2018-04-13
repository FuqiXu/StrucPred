cd BioinfoProject/KB8025/scripts/SequencesToblast
for i in *.fasta
do
	echo "PSI-BLAST $i"
	psiblast -db swissprot/uniprot_sprot.fasta -query $i -outfmt 10 -num_iterations 3 -num_threads 8 -comp_based_stats 1 -out_ascii_pssm result/$i.out
	echo "finish $i"
done

echo "Completed"
