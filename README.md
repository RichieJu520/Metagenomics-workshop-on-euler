# Metagenomics-workshop-on-euler
A workshop for microbial metagenomics and bioinformatics practice on the euler server of ETH Zurich


91_otus.fasta.gz         --- 91%-identity OTUs of 16S rRNA gene Greengenes Database  http://greengenes.secondgenome.com/downloads
SARG_20170328.fasta.gz   --- Structure antibiotic resistance gene database (ARDB) downloaded from https://smile.hku.hk/SARGs

D10.R1_1.fq.gz           --- shot-gun metagenomes with paired end reads 1 for sample D10.R1
D10.R1_2.fq.gz           --- shot-gun metagenomes with paired end reads 2 for sample D10.R1
D10.R2_1.fq.gz           --- shot-gun metagenomes with paired end reads 1 for sample D10.R2
D10.R2_2.fq.gz           --- shot-gun metagenomes with paired end reads 2 for sample D10.R2

N50_GC.py                --- a python script for calculating N50, average, minimum, and maximum length of sequences in multiple fastas
submit.blastp.jobarrays.cmds.lsf 
                         --- a lsf file for bsub submission of a blast jobarrays to server
calc.coverage.in.bam.depth.pl  (copyright: Mads Albertsen)
                         --- a perl script for calculating coverage from the depth profile of a sorted bam file (sorted by samtools)               
