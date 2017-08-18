# pu-variants
Finds variant sites from samtools mpileup output
To make:
>make pu-variants
To run:
>samtools mpileup -s [other filtering options, as necessary] -f [fasta file] input.bam | ./pu-variants -l [low cov] -h [high cov] -H
To set required number of observations of alleles on both strands, set MIN_ALLELE_BOTH_STRANDS in pu-variant.c source and recompile
