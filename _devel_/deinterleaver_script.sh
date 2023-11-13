for INF in *.fasta; do
cat $INF | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | tail -n +2 > ${INF%.fasta*}_deint.fasta; 
done
