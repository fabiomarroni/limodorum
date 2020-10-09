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


#Download and gunzip the Uniprot GO:
cd ${PROJDIR}/GO
wget -c ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
gunzip goa_uniprot_all.gaf.gz

cd ${FUNCDIR}/scripts/logs
BLASTDB=/iga/biodb/ncbi/blastdb/latest/swissprot
MYQ=${PROJDIR}/full_trinity_assembly/Trinity_new.fasta
MYB=${PROJDIR}/GO/Trinity_blast.out
MYGO=${PROJDIR}/GO/goa_uniprot_all.gaf
blast1=`echo "module load aligners/blast/latest; blastx -db $BLASTDB -query $MYQ -num_threads 12\
 -max_target_seqs 2 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -evalue 1e-10 -out $MYB" | \
qsub -N blast1 -l vmem=24G,walltime=168:00:00,nodes=1:ppn=12`


#Run on the second part of file after the first was terminated by walltime 
MYQ2=${PROJDIR}/full_trinity_assembly/Trinity_new_2.fasta
MYB2=${PROJDIR}/GO/Trinity_blast_2.out
blast1=`echo "module load aligners/blast/latest; blastx -db $BLASTDB -query $MYQ2 -num_threads 12\
 -max_target_seqs 2 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -evalue 1e-10 -out $MYB2" | \
qsub -N blast1 -l vmem=24G,walltime=168:00:00,nodes=1:ppn=12`

FULLBLAST=${PROJDIR}/GO/Trinity_full_blast.out
#Merge the two files
echo "cat $MYB $MYB2 > $FULLBLAST" | qsub -N merge -W depend=afterok:$blast1



#Here I start to prepare scripts for parsing the blastx output and searching for the Uniprot/Swissprot IDs in the uniprot GO file.
##########################################################################
#Extract the Swissprot ID from blast output and write to file
##########################################################################
BSWISS=${PROJDIR}/GO/blasted_swissprot_ID.txt
cut -f2 $FULLBLAST | cut -d"|" -f4| cut -d"." -f1 | sort | uniq > $BSWISS

#Given a swissprot ID search it in the GO file
#Test command, no need to run it
#awk '$2=="A0A000"' $MYGO



##########################################################################
#Read the gene_id in $BSWISS and return them if they have a correspondence in the second column of $MYGO
##########################################################################

MYMAPPING=${PROJDIR}/GO/transcripts_GO.txt
#This might be very slow, and we might rather run the command below on the split files.
awk 'FNR==NR{k[$1]=1;next;} k[$2]' $BSWISS $MYGO > $MYMAPPING 



#Subset the large go file. Note the exit command that avoids reading all the lines after 100000
# 608472856
GOSUBDIR=${PROJDIR}/GO/split_GO
mkdir -p $GOSUBDIR
MYLINE=(0 100000000 200000000 300000000 400000000 500000000 600000000 608472856)
for((i=1;i<${#MYLINE[@]};i++))
do
	START=$[MYLINE[$i-1]+1]
	END=$[MYLINE[$i]]   
	echo $START
	echo $END	
	echo "print from $[i-1]: $START to $i: $END"
	OUTPREF=$(basename $MYGO)	
	awk -v mystart="${START}" -v myend="${END}" 'BEGIN {print myend}'
	awk -v mystart="${START}" -v myend="${END}" 'NR > myend { exit } NR >= mystart && NR <= myend' $MYGO > $GOSUBDIR/${OUTPREF/.gaf/_${i}.gaf}
done

#Loop over the "small" GO files and match the genes identified by blastx results against them.
#We will obtain the same number of files as the subset, each one reporting the match of the full set of transcript blastx
#against the subset of GO. We will then have to merge them, probably reformat them and then merge them with the file of the DE analysis.
MATCHDIR=${PROJDIR}/GO/matching_GO
mkdir -p $MATCHDIR
for SUBGO in $GOSUBDIR/*gaf
#for SUBGO in $GOSUBDIR/*7.gaf
do
echo "Reading from file $SUBGO"
SMATCH=$(basename $SUBGO)
SMATCH=${SMATCH/.gaf/.txt}
echo "Writing to file ${MATCHDIR}/${SMATCH}"
#awk 'FNR==NR{k[$1]=1;next;} FNR==k[$2]' $BSWISS $SUBGO
awk 'FNR==NR{k[$1]=1;next;} k[$2]' $BSWISS $SUBGO > ${MATCHDIR}/${SMATCH}
#grep -f $BSWISS $SUBGO > ${MATCHDIR}/${SMATCH}
done

#awk 'FNR==NR{k[$1]=1;next;} k[$2]' $BSWISS test.gaf

#Reformat GO file, so that each gene has only one row and all GO terms are separated by ";"
FORMGO=${PROJDIR}/GO/formatted_transcripts_GO.txt
Rscript ${FUNCDIR}/functions/limodorum/03_reformat_GO.r \
-I $MYMAPPING \
-O $FORMGO


#################################
#OLD STUFF BELOW (This  was done by testing GO and KEGG on the full dataset just fo rfun, but we then preferred to perform analysis separated for plants and fungi)
#################################

#Merge GO terms (in $MYMAPPING) to the DE file and produce a summary table aggregating data by phylum
RESP=${PROJDIR}/DE/DE_clean_phyla.txt
OUTFILE=${PROJDIR}/DE/DE_GO.txt
SUMMARY=${PROJDIR}/Tables/phylum_summary.txt
Rscript ${FUNCDIR}/functions/limodorum/04_merge_blast_GO_DE.r \
-G $FORMGO \
-B $FULLBLAST \
-D $RESP \
-S $SUMMARY \
-O $OUTFILE

#ADD KEGG TERMS
INFILE=${PROJDIR}/DE/DE_GO.txt
OUTFILE=${PROJDIR}/DE/DE_GO_KEGG.txt
Rscript ${FUNCDIR}/functions/limodorum/05_assign_kegg.r \
-I $INFILE \
-O $OUTFILE

FUNGIFILE=${PROJDIR}/DE/DE_clean_ABmycota.txt
PLANTFILE=${PROJDIR}/DE/DE_clean_Streptophyta.txt
#Create DE files including only Basidiomycota and Streptophyta
awk -F"\t" 'NR == 1 || $11=="Basidiomycota"||$11=="Ascomycota" {print $0}' $OUTFILE > $FUNGIFILE
awk -F"\t" 'NR == 1 || $11=="Streptophyta" {print $0}' $OUTFILE > $PLANTFILE





#Run GO enrichment test using topGO

#ON Fungi
for ONTO in BP CC MF
do
OUTFILE=${PROJDIR}/DE/DE_ABmycota_GO_enrich_${ONTO}.txt
Rscript ${FUNCDIR}/functions/limodorum/06_enrichment.r -G $FORMGO -I $FUNGIFILE -N $ONTO -O $OUTFILE
done

#ON Plants
for ONTO in BP CC MF
do
OUTFILE=${PROJDIR}/DE/DE_ABmycota_GO_enrich_${ONTO}.txt
Rscript ${FUNCDIR}/functions/limodorum/06_enrichment.r -G $FORMGO -I $PLANTFILE -N $ONTO -O $OUTFILE
done





