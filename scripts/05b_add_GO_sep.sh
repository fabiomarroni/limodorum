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

#This is similar to 05, but we use DE files obtained by keeping fungi and plants separate. This step (adding phylum) is not necessary since we already know who is who, 
# but is needed to properly reformat DE files!

#I only deleted the once only tasks, such as downloading databases or performing blastx


FULLBLAST=${PROJDIR}/GO/Trinity_full_blast.out
BSWISS=${PROJDIR}/GO/blasted_swissprot_ID.txt
MYMAPPING=${PROJDIR}/GO/transcripts_GO.txt
FORMGO=${PROJDIR}/GO/formatted_transcripts_GO.txt


for ORG in AB plants
do
#Again the trick for bringing with us the average expression for LM and LS based on raw or normalized counts
for count in norm raw
do
#Merge GO terms (in $MYMAPPING) to the DE file and produce a summary table aggregating data by phylum
RESP=${PROJDIR}/DE_${ORG}/${ORG}_DE_${count}_2020_new.txt
OUTFILE=${PROJDIR}/DE_${ORG}/${ORG}_DE_${count}_GO_2020.txt
Rscript ${FUNCDIR}/functions/limodorum/04c_merge_blast_GO_DE_no_phylum.r \
-G $FORMGO \
-B $FULLBLAST \
-D $RESP \
-O $OUTFILE

#ADD KEGG TERMS
INFILE=${PROJDIR}/DE_${ORG}/${ORG}_DE_${count}_GO_2020.txt
KEGGDE=${PROJDIR}/DE_${ORG}/${ORG}_DE_${count}_GO_KEGG_2020.txt
Rscript ${FUNCDIR}/functions/limodorum/05_assign_kegg.r \
-I $INFILE \
-O $KEGGDE
done
done


#Run GO enrichment test using topGO (output is writeen as .txt and as .xlsx)
for ORG in AB plants
do
KEGGDE=${PROJDIR}/DE_${ORG}/DE_GO_KEGG.txt
for ONTO in BP CC MF
    do
    OUTFILE=${PROJDIR}/DE_${ORG}/DE_GO_enrich_${ONTO}_${ORG}.txt
    Rscript ${FUNCDIR}/functions/limodorum/07_enrichment_sep.r -G $FORMGO -I $KEGGDE -N $ONTO -O $OUTFILE
    done
done


#Run GO enrichment WITHOUT using topGO
for ORG in AB plants
do
KEGGDE=${PROJDIR}/DE_${ORG}/DE_GO_KEGG.txt
for ONTO in BP CC MF
    do
    OUTFILE=${PROJDIR}/DE_${ORG}/DE_GO_notopGO_enrich_${ONTO}_${ORG}.txt
    Rscript ${FUNCDIR}/functions/limodorum/07a_enrich_GO_sep.r -I $KEGGDE -N $ONTO -O $OUTFILE
    done
done



#FOR GENES WITH KEGG, ASSIGN PATHWAY (To prevent possible format conflicts, it is better to run this after we tested GO enrichment)
for count in norm raw
do
for ORG in AB plants
do
    INFILE=${PROJDIR}/DE_${ORG}/${ORG}_DE_${count}_GO_KEGG_2020.txt
    OUTFILE=${PROJDIR}/DE_${ORG}/${ORG}_DE_${count}_GO_KEGG_path_2020.txt
    Rscript ${FUNCDIR}/functions/limodorum/08_kegg_term_2_path.r \
    -I $INFILE \
    -O $OUTFILE
done
done

#If an old GO descritpion column is left, remove it (since it may cause heart attack seeing Homo sapiens in fungi!!)
#Also if LM and LS columns are left remove them, because they are WRONG!!!! We will add them at the end!
for count in norm raw
do
for ORG in AB plants
do
    OUTFILE=${PROJDIR}/DE_${ORG}/${ORG}_DE_${count}_GO_KEGG_path_2020.txt
    NODESCFILE=${PROJDIR}/DE_${ORG}/${ORG}_DE_${count}_GO_KEGG_clean_2020.txt
    Rscript ${FUNCDIR}/functions/limodorum/16_remove_descr.r \
    -I $OUTFILE \
	-R description,LM,LS \
    -O $NODESCFILE
done
done



#TEST ENRICHMENT IN KEGG PATHWAYS
count=raw
for ORG in AB plants
do
    INFILE=${PROJDIR}/DE_${ORG}/${ORG}_DE_${count}_GO_KEGG_clean_2020.txt
    OUTFILE=${PROJDIR}/DE_${ORG}/KEGG_enrich_${ORG}.txt
    Rscript ${FUNCDIR}/functions/limodorum/09_kegg_enrich.r \
    -I $INFILE \
    -O $OUTFILE
done



#Plot GO enrichment
Rscript ${FUNCDIR}/functions/limodorum/12_plot_GO_enrich.r \
-F ${PROJDIR}/DE_AB/ \
-P ${PROJDIR}/DE_plants/ \
-O ${PROJDIR}/plots/GO_enrich.png

#Plot GO enrichment NOT using topGO results (we use a fisher test as we did in KEGG)
Rscript ${FUNCDIR}/functions/limodorum/12a_plot_GO_notopGO_enrich.r
-F ${PROJDIR}/DE_AB/ \
-P ${PROJDIR}/DE_plants/ \
-O ${PROJDIR}/plots/GO_notopGO_enrich.png \
-T ${PROJDIR}/Tables

#Plot KEGG enrichment
Rscript ${FUNCDIR}/functions/limodorum/13_plot_KEGG_enrich.r \
-F ${PROJDIR}/DE_AB/KEGG_enrich_AB.txt \
-P ${PROJDIR}/DE_plants/KEGG_enrich_plants.txt \
-O ${PROJDIR}/plots/KEGG_enrich_2020.png




