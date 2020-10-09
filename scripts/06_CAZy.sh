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

#Download the CAZy database fromd DBcan
CAZYDIR=${PROJDIR}/CAZy
mkdir -p $CAZYDIR
cd $CAZYDIR
wget -c http://bcb.unl.edu/dbCAN2/download/CAZyDB.07312019.fa

# Then I realized that I preferred using the DBcan annotation tool
# We need to upload the file, which must be at most 20MB in size
cd $CAZYDIR/ref_dbcan
#Split with pyfasta
#We can also split using awk, but I didn't try that... Might like to play with that! (https://www.biostars.org/p/13270/)

pyfasta split -n 30 ../../Reference/Trinity.fasta
for aaa in ../../Reference/Trinity.??.fasta ; do echo $aaa; mv $aaa .; done
 #I then upload the files to http://bcb.unl.edu/dbCAN2/blast.php selecting all the analysis methods

#Now we merge the results and provide a report of the cazy results
for count in norm raw
do
for aaa in AB plants
do
Rscript ${FUNCDIR}/functions/limodorum/14_add_CAZY.r \
-I ${PROJDIR}/CAZy/dbcan2 \
-D ${PROJDIR}/DE_${aaa}/${aaa}_DE_${count}_GO_KEGG_clean_2020.txt \
-O ${PROJDIR}/DE_${aaa}/${aaa}_DE_${count}_GO_KEGG_path_cazy.txt
done
done

#Add Reads number (overwrite file)
#This is a horrible patch I made a posteriori to add read numbers for LM and LS to the annotated files.

for count in norm raw
do
for aaa in AB plants
do
Rscript ${FUNCDIR}/functions/limodorum/17_add_reads_counts.r \
-I ${PROJDIR}/DE_${aaa}/${aaa}_DE_${count}_GO_KEGG_path_cazy.txt \
-D ${PROJDIR}/DE_${aaa}/${aaa}_DE_${count}_2020_new.txt \
-O ${PROJDIR}/DE_${aaa}/${aaa}_DE_${count}_GO_KEGG_path_cazy.txt
done
done




#Perform enrichment analysis and plots of CAZY terms
count=raw
Rscript ${FUNCDIR}/functions/limodorum/15_enrich_CAZY.r \
-P ${PROJDIR}/DE_plants/plants_DE_${count}_GO_KEGG_path_cazy.txt \
-F ${PROJDIR}/DE_AB/AB_DE_${count}_GO_KEGG_path_cazy.txt \
-O ${PROJDIR}/plots/DE_GO_KEGG_path_cazy_enrich_2020.png
