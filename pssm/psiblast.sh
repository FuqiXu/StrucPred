for i in Sequences/*.fasta
do
	echo "PSI-BLAST $i"
	psiblast -db uniprot_sprot.fasta -query $i -outfmt 10 -num_iterations 3 -num_threads 8 -comp_based_stats 1 -out_ascii_pssm psipssm/$i.pssm
	echo "finish $i"
done

echo "Completed"
