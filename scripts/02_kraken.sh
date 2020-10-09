#!/bin/bash

##############################################
#Variables you will need to change
###############################################
PROJDIR=/path/to/results  #Everything will be saved as a subfolder of this
SOFTDIR=/path/to/software #Assuming you have a common repository for software (only needed if the software packages I used are not in your path)
FUNCDIR=/path/to/scripts/and/functions #I suggest you create inside this folder a folder called "functions/limodorum" in which you will save all the R functions provided.
										# You also need to create a folder scripts/limodorum in which you need to put al the shell and python scripts
ANACDIR=/path/to/anaconda/bins #I needed it to quickly access my conda environments
KRAKBIN=/path/to/kraken/executables # Only needed if kraken is not in your path nor in SOFTDIR (I am a mess, and it was in its own folder!)
KDB=/path/to/kraken/db #No comments
KSCRIPT=/path/to/kraken/helper/scripts #
TAXKIT=/path/to/taxonkit/bin        #Only needed if taxonkit is not in your repository specified with SOFTDIR
##############################################
#End of variables you will need to change
###############################################


INPUT=${PROJDIR}/full_trinity_assembly/Trinity.fasta
KRES=${PROJDIR}/kraken_nt
outpref=Trinity
# PRE Execution
#s1 - classify
mkdir -p ${KRES}/logs
krak1=`echo "cd ${KRES}/logs; export TMPDIR=${KRES}; 
module use --append /iga/scripts/dev_modules/modules; \
${KRAKBIN}/kraken2 --threads 16 --db $KDB $INPUT --unclassified-out ${KRES}/unclassified.fasta --output ${KRES}/${outpref}.kraken --use-names --report ${KRES}/${outpref}.kraken.report.txt >${KRES}/logs/krak.out 2>${KRES}/logs/krak.err"| qsub -N krak1_${outpref} -l vmem=220G,walltime=168:00:00,nodes=1:ppn=16`
#s2 - reports
echo "cd ${KRES}; \
export TMPDIR=${KRES}; \
module use --append /iga/scripts/dev_modules/modules; \
module load dev/krona/2.6; \
python  ${KSCRIPT}/kraken2txt.py ${KRES}/${outpref}.kraken.report.txt ${KRES}/${outpref}.kraken.krona_table.txt; \
ktImportText ${KRES}/${outpref}.kraken.krona_table.txt,${outpref} -o ${outpref}.kraken_single-chart.html >${KRES}/logs/map.out 2>${KRES}/logs/map.err" | qsub -N krak2_${outpref} -l vmem=32G,walltime=24:00:00,nodes=1:ppn=16 -W depend=afterok:$krak1



#Run KRAKEN on Reads, mostly to check to which fungus we are mostly assigning reads.
READ_LEN=150

cd ${FUNCDIR}/scripts/logs
PRES=${PROJDIR}
KRES=${PROJDIR}/kraken_reads
mkdir -p $KRES/logs
for SAMPLE in LM3 LM4 LM6 LS2 LS3 LS6
do
trimmed_dir=${PRES}/Sample_${SAMPLE}/merge_trimmed 
READ1=${trimmed_dir}/${SAMPLE}_R1.fastq.gz
READ2=${trimmed_dir}/${SAMPLE}_R2.fastq.gz
READ2=${READ1/_norRNA_1/_norRNA_2}
pref=$SAMPLE
step1=`echo "cd ${KRES}; export TMPDIR=${KRES}; module use --append /iga/scripts/dev_modules/modules; \
${KRAKBIN}/kraken2 --threads 16 --paired --gzip-compressed --db $KDB ${READ1} ${READ2} --output ${KRES}/${pref}.kraken --use-names --report ${KRES}/${pref}.kraken.report.txt >${KRES}/logs/map.out 2>${KRES}/logs/map.err"| qsub -N s1_${pref} -l vmem=150G,walltime=168:00:00,nodes=1:ppn=16`
#s2 - reports
echo "cd ${KRES}; \
export TMPDIR=${KRES}; \
module use --append /iga/scripts/dev_modules/modules; \
module load dev/krona/2.6; \
python ${KSCRIPT}/kraken2txt.py ${KRES}/${pref}.kraken.report.txt ${KRES}/${pref}.kraken.krona_table.txt; \
ktImportText ${KRES}/${pref}.kraken.krona_table.txt,${pref} -o ${pref}.kraken_single-chart.html >${KRES}/logs/map.out 2>${KRES}/logs/map.err" | qsub -N s2_${pref} -l vmem=32G,walltime=24:00:00,nodes=1:ppn=16 -W depend=afterok:$step1
done

#Use blastx to classify transcripts, just a second check
#!/bin/bash
module load aligners/blast/latest
INPUT=${PROJDIR}/full_trinity_assembly/Trinity.fasta
OUTPUT=${PROJDIR}/full_trinity_assembly/Trinity_blast_nr.out
NT=16
echo "module load aligners/blast/latest; blastx -query $INPUT -db /iga/biodb/ncbi/blastdb/latest/nr -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -evalue 1e-10 -num_alignments 1 -num_threads $NT -out $OUTPUT" | qsub -N balnr -l vmem=32G,walltime=168:00:00,nodes=1:ppn=$NT

module load aligners/blast/latest
INPUT=${PROJDIR}/full_trinity_assembly/Trinity.fasta
OUTPUT=${PROJDIR}/full_trinity_assembly/Trinity_blast_nt.out
echo "module load aligners/blast/latest;blastn -query $INPUT -db /iga/biodb/ncbi/blastdb/latest/nt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -evalue 1e-10 -num_alignments 1 -num_threads $NT -out $OUTPUT " | qsub -N balnr -l vmem=32G,walltime=168:00:00,nodes=1:ppn=$NT

