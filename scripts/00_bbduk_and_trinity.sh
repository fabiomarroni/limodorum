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



#Run bbduk on all samples
for MYDIR in ${PROJDIR}/Sample_*
do
for READ1 in ${MYDIR}/*_R1_*
do
READ1=$(basename $READ1)
READ2=${READ1/_R1_/_R2_}
echo $READ1
echo $READ2
TRIMDIR=${MYDIR}/trimmed
${SOFTDIR}/bbmap/bbduk.sh in1=${MYDIR}/${READ1} in2=${MYDIR}/${READ2} out1=${TRIMDIR}/${READ1} out2=${TRIMDIR}/${READ2} ref=${SOFTDIR}/bbmap/resources/adapters.fa ktrim=r ktrim=l qtrim=rl trimq=10 tpe tbo
done
done


#On full set, trimmed reads
LOGDIR=${FUNCDIR}/scripts/logs
NCORES=16
PROJDIR=${PROJDIR}/
LM3=${PROJDIR}/Sample_LM3/trimmed
LM4=${PROJDIR}/Sample_LM4/trimmed
LM6=${PROJDIR}/Sample_LM6/trimmed
LS2=${PROJDIR}/Sample_LS2/trimmed
LS3=${PROJDIR}/Sample_LS3/trimmed
LS6=${PROJDIR}/Sample_LS6/trimmed
LREADS=${LM3}/LM3_ATGTCA_L002_R1_001.fastq.gz,${LM3}/LM3_ATGTCA_L002_R1_002.fastq.gz,${LM3}/LM3_ATGTCA_L002_R1_003.fastq.gz,${LM3}/LM3_ATGTCA_L002_R1_004.fastq.gz,${LM4}/LM4_CCGTCC_L002_R1_001.fastq.gz,${LM4}/LM4_CCGTCC_L002_R1_002.fastq.gz,${LM4}/LM4_CCGTCC_L002_R1_003.fastq.gz,${LM4}/LM4_CCGTCC_L002_R1_004.fastq.gz,${LM6}/LM6_GTCCGC_L002_R1_001.fastq.gz,${LM6}/LM6_GTCCGC_L002_R1_002.fastq.gz,${LM6}/LM6_GTCCGC_L002_R1_003.fastq.gz,${LM6}/LM6_GTCCGC_L002_R1_004.fastq.gz,${LS2}/LS2_AGTCAA_L002_R1_001.fastq.gz,${LS2}/LS2_AGTCAA_L002_R1_002.fastq.gz,${LS2}/LS2_AGTCAA_L002_R1_003.fastq.gz,${LS2}/LS2_AGTCAA_L002_R1_004.fastq.gz,${LS3}/LS3_AGTTCC_L002_R1_001.fastq.gz,${LS3}/LS3_AGTTCC_L002_R1_002.fastq.gz,${LS3}/LS3_AGTTCC_L002_R1_003.fastq.gz,${LS3}/LS3_AGTTCC_L002_R1_004.fastq.gz,${LS6}/LS6_GTGAAA_L002_R1_001.fastq.gz,${LS6}/LS6_GTGAAA_L002_R1_002.fastq.gz,${LS6}/LS6_GTGAAA_L002_R1_003.fastq.gz,${LS6}/LS6_GTGAAA_L002_R1_004.fastq.gz
RREADS=$(echo $LREADS | sed -e 's/_R1_/_R2_/'g)

#echo ${READNAMES}                 
echo "export PATH=$PATH:${ANACDIR}; source activate trinity_2.5.1; Trinity --seqType fq --max_memory 120G --left ${LREADS} --right ${RREADS} --output ${PROJDIR}/full_trinity_assembly --CPU $NCORES >${LOGDIR}/trinity.out 2>${LOGDIR}/trinity.err" | qsub -N trinity -l vmem=128G,walltime=168:00:00,nodes=1:ppn=$NCORES


#Try to run trinity assuming stranded reads
#Also remember that also hisat2, stringtie and htseq count probably need to be run using strand information


echo "export PATH=$PATH:${ANACDIR}; source activate trinity_2.5.1; Trinity --seqType fq --max_memory 120G --left ${LREADS} --right ${RREADS} --output ${PROJDIR}/full_trinity_assembly_RF --CPU $NCORES --SS_lib_type RF >${LOGDIR}/trinity.out 2>${LOGDIR}/trinity.err" | qsub -N trinity -l vmem=128G,walltime=168:00:00,nodes=1:ppn=$NCORES
echo "export PATH=$PATH:${ANACDIR}; source activate trinity_2.5.1; Trinity --seqType fq --max_memory 120G --left ${LREADS} --right ${RREADS} --output ${PROJDIR}/full_trinity_assembly_FR --CPU $NCORES --SS_lib_type FR >${LOGDIR}/trinity.out 2>${LOGDIR}/trinity.err" | qsub -N trinity -l vmem=128G,walltime=168:00:00,nodes=1:ppn=$NCORES



#Merge reads for each sample (I could have done that BEFORE running trinity, but I forgot) to subsequently
#perform differential expression analysis.

for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
SAMPLEINDIR=${PROJDIR}/Sample_${SAMPLENAME}/trimmed
SAMPLEOUTDIR=${PROJDIR}/Sample_${SAMPLENAME}/merge_trimmed
mkdir -p $SAMPLEOUTDIR
#Put file names in an array
declare -a READS1=$(ls $SAMPLEINDIR/*_R1_*fastq.gz)
declare -a READS2=$(ls $SAMPLEINDIR/*_R2_*fastq.gz)
#Cat all files listed in the array to an output file
cat ${READS1[@]} > ${SAMPLEOUTDIR}/${SAMPLENAME}_R1.fastq.gz
cat ${READS2[@]} > ${SAMPLEOUTDIR}/${SAMPLENAME}_R2.fastq.gz
done


#Count raw reads
MYARRAY=()
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
SAMPLEOUTDIR=${PROJDIR}/Sample_${SAMPLENAME}
READS1=$(ls $SAMPLEOUTDIR/*_R1_*.fastq.gz)
MYLINES=$(zcat $READS1 | wc -l)
MYOUT=$( echo $READS1 $(($MYLINES / 4 )))
MYARRAY+=("$MYOUT")
echo $(($MYLINES / 4 ))
done
printf "%s\n" "${MYARRAY[@]}" > ${PROJDIR}/raw_reads.txt


#Count merged trimmed reads
MYARRAY=()
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
SAMPLEOUTDIR=${PROJDIR}/Sample_${SAMPLENAME}/merge_trimmed
READS1=$(ls $SAMPLEOUTDIR/*_R1.fastq.gz)
MYLINES=$(zcat $READS1 | wc -l)
MYOUT=$( echo $READS1 $(($MYLINES / 4 )))
MYARRAY+=("$MYOUT")
echo $(($MYLINES / 4 ))
done
printf "%s\n" "${MYARRAY[@]}" > ${PROJDIR}/merged_trimmed_reads.txt


