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


#Extract ids of streptophyta, basidiomycota and ascomycota
${TAXKIT}/taxonkit list --ids 35493 --indent "" --show-name --data-dir /iga/biodb/metagenomics/krakenuniq/taxonomy \
> ${PROJDIR}/DE/Streptophyta_names.txt
${TAXKIT}/taxonkit list --ids 35493 --indent "" --data-dir /iga/biodb/metagenomics/krakenuniq/taxonomy \
> ${PROJDIR}/DE/Streptophyta.txt
${TAXKIT}/taxonkit list --ids 4890,5204 --indent "" --show-name --data-dir /iga/biodb/metagenomics/krakenuniq/taxonomy \
> ${PROJDIR}/DE/AB_names.txt
${TAXKIT}/taxonkit list --ids 4890,5204 --indent "" --data-dir /iga/biodb/metagenomics/krakenuniq/taxonomy \
> ${PROJDIR}/DE/AB.txt


#Get statistics of unique alignments in plants and AB
for ORG in AB plants
do
declare -a MYARRAY=()
COUNTDIR=${PROJDIR}/htseq_${ORG}
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
MYCOUNT=$(awk ' {sum +=$2} END {print sum}' $COUNTDIR/${SAMPLENAME}.txt)
MYARRAY+=("$SAMPLENAME $MYCOUNT")
done
printf "%s\n" "${MYARRAY[@]}" > ${COUNTDIR}/total.txt
done

#Get statistics of unique alignments in the full dataset
declare -a MYARRAY=()
COUNTDIR=${PROJDIR}/htseq
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
MYCOUNT=$(grep -v "__" $COUNTDIR/${SAMPLENAME}.txt | awk ' {sum +=$2} END {print sum}' )
MYARRAY+=("$SAMPLENAME $MYCOUNT")
done
printf "%s\n" "${MYARRAY[@]}" > ${COUNTDIR}/total.txt




#Read DE files with assignments to species and using the taxonid output separate streptophyta from fungi and discard all others
Rscript ${FUNCDIR}/functions/limodorum/02_sep_plants_AB.r \
-D ${PROJDIR}/DE/DE_new.txt \
-T ${PROJDIR}/DE/Streptophyta.txt \
-O ${PROJDIR}/DE_plants/DE_new.txt

Rscript ${FUNCDIR}/functions/limodorum/02_sep_plants_AB.r \
-D ${PROJDIR}/DE/DE_new.txt \
-T ${PROJDIR}/DE/AB.txt \
-O ${PROJDIR}/DE_AB/DE_new.txt


#I copy the file from streptophyta to plants becaunse plants is the common name I am going to use
cp ${PROJDIR}/DE/Streptophyta.txt ${PROJDIR}/DE/plants.txt
#Plot the number of reads attributed by kraken to each fungal or plant genus.
for ORG in AB plants
do
Rscript ${FUNCDIR}/functions/limodorum/10_plot_species_kraken.r \
-I ${PROJDIR}/kraken_reads \
-T ${PROJDIR}/DE/${ORG}.txt \
-L G \
-O ${PROJDIR}/DE_${ORG}/${ORG}_genus.png
done

#Plot fungal genera separately for LM and LS
Rscript ${FUNCDIR}/functions/limodorum/10a_plot_fungi_kraken.r \
-I ${PROJDIR}/kraken_reads \
-T ${PROJDIR}/DE/AB.txt \
-L G \
-S FALSE \
-O ${PROJDIR}/DE_AB/Fig_4_AB_genus.png


#Use blast to further classify Fungi. 
#We need to check what's happening
F_TRANSCRIPT=${PROJDIR}/DE_AB/AB_sequence_names.txt
cut -f8 ${PROJDIR}/DE_AB/DE_new.txt | sed 1d > $F_TRANSCRIPT
MYQ=${PROJDIR}/full_trinity_assembly/Trinity.fasta
F_SEQ=${PROJDIR}/full_trinity_assembly/AB.fasta
#Extract fungal transcripts
module load sw/bio/seqtk/default
seqtk subseq $MYQ $F_TRANSCRIPT > $F_SEQ
#Run blastn
BLASTDB=/iga/biodb/ncbi/blastdb/latest/nt
MYB2=${PROJDIR}/DE_AB/AB_blast.out
CPU=24
#I also write taxon id and species, so that I don't have to look them up later
blast1=`echo "module load aligners/blast/latest; blastn -db $BLASTDB -query $F_SEQ -num_threads $CPU\
 -max_target_seqs 2 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames' -evalue 1e-10 -out $MYB2" | \
qsub -N blast1 -l vmem=24G,walltime=168:00:00,nodes=1:ppn=$CPU`


#Use blast to map against SILVA
module load aligners/blast/latest
SILVADB=${PROJDIR}/Reference/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta
makeblastdb -in $SILVADB -dbtype nucl
CPU=8
MYSIL=${PROJDIR}/DE_AB/AB_silva.out
blast1=`echo "module load aligners/blast/latest; blastn -db $SILVADB -query $F_SEQ -num_threads $CPU\
 -max_target_seqs 2 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames' -evalue 1e-10 -out $MYSIL" | \
qsub -N blast1 -l vmem=24G,walltime=168:00:00,nodes=1:ppn=$CPU`

