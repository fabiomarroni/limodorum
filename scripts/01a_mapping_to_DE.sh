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


#Here, several attempts were made, to check if paired reads were stranded or not because no information was available.
#Maybe I could have found a quicker check, but I did a lot of check and at the end
#I found out that the reads were not stranded (and this corresponded to what my collaborators thought was the case, so we were happy!)
#I leave everything just in case we need to run some more test iwht strandedness.

#If your software doesn't runt hrough modules you don't need this.
#If you module runs thorugh modules you might need to load the correct module
module load sw/aligners/hisat2/2.0.4
module load sw/bio/samtools/0.1.18

cp ${PROJDIR}/full_trinity_assembly/Trinity.fasta ${PROJDIR}/Reference

cp ${PROJDIR}/full_trinity_assembly_FR/Trinity.fasta ${PROJDIR}/Reference/Trinity_FR.fasta
cp ${PROJDIR}/full_trinity_assembly_RF/Trinity.fasta ${PROJDIR}/Reference/Trinity_RF.fasta


echo "module load sw/aligners/hisat2/2.0.4; hisat2-build ${PROJDIR}/Reference/Trinity.fasta ${PROJDIR}/Reference/Trinity" | qsub -N his_build -l vmem=32G,walltime=10:00:00,nodes=1:ppn=4
#Running on Oct 9th 2019
echo "module load sw/aligners/hisat2/2.0.4; hisat2-build ${PROJDIR}/Reference/Trinity_FR.fasta ${PROJDIR}/Reference/Trinity_FR" | qsub -N his_build -l vmem=32G,walltime=10:00:00,nodes=1:ppn=4
echo "module load sw/aligners/hisat2/2.0.4; hisat2-build ${PROJDIR}/Reference/Trinity_RF.fasta ${PROJDIR}/Reference/Trinity_RF" | qsub -N his_build -l vmem=32G,walltime=10:00:00,nodes=1:ppn=4



#Run hisat2 in all samples assuming unstranded libraries
OUTDIR=${PROJDIR}/h2_align
mkdir -p $OUTDIR
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
MY_IN1=${PROJDIR}/Sample_${SAMPLENAME}/merge_trimmed/*R1.fastq.gz
MY_IN2=${PROJDIR}/Sample_${SAMPLENAME}/merge_trimmed/*R2.fastq.gz
echo "module load sw/aligners/hisat2/2.0.4;hisat2 ${PROJDIR}/Reference/Trinity --no-unal -1 $MY_IN1 -2 $MY_IN2 -S ${OUTDIR}/${SAMPLENAME}.sam" | qsub -N h2_align -l vmem=32G,walltime=24:00:00,nodes=1:ppn=4 
done
#Run hisat2 in all samples assuming stranded libraries (we try both)
UPDIR=${PROJDIR}/h2_align
NCORES=2
for STRAND in RF FR
do
OUTDIR=${UPDIR}_${STRAND}
echo $OUTDIR
mkdir -p $OUTDIR
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
MY_IN1=${PROJDIR}/Sample_${SAMPLENAME}/merge_trimmed/*R1.fastq.gz
MY_IN2=${PROJDIR}/Sample_${SAMPLENAME}/merge_trimmed/*R2.fastq.gz
echo "module load sw/aligners/hisat2/2.0.4;hisat2 ${PROJDIR}/Reference/Trinity_${STRAND} --no-unal -1 $MY_IN1 -2 $MY_IN2 --${STRAND,,} -p $NCORES -S ${OUTDIR}/${SAMPLENAME}.sam" | qsub -N h2_align -l vmem=32G,walltime=24:00:00,nodes=1:ppn=$NCORES
done
done

#Convert alignments to bam assuming unstrandedness
module load sw/bio/samtools/0.1.18
MYARRAY=()
SAMPLEOUTDIR=${PROJDIR}/h2_align
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
samtools view -bS ${SAMPLEOUTDIR}/${SAMPLENAME}.sam -o ${SAMPLEOUTDIR}/${SAMPLENAME}.bam 
rm ${SAMPLEOUTDIR}/${SAMPLENAME}.sam
done

#Convert alignments to bam assuming strandedness (we try both)
module load sw/bio/samtools/0.1.18
for STRAND in RF FR
do
SAMPLEOUTDIR=${PROJDIR}/h2_align_${STRAND}
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
samtools view -bS ${SAMPLEOUTDIR}/${SAMPLENAME}.sam -o ${SAMPLEOUTDIR}/${SAMPLENAME}.bam 
rm ${SAMPLEOUTDIR}/${SAMPLENAME}.sam
done
done


#READ COUNTS. Some of them will be put in Table 1 of the manuscript

#Count uniquely aligned reads (uniquely mapping reads do not have ZS:i field) assuming no strandedness
#To count the number of uniquely mapped reads 
module load sw/bio/samtools/0.1.18
MYARRAY=()
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
SAMPLEOUTDIR=${PROJDIR}/h2_align
MYOUT=$( samtools view ${SAMPLEOUTDIR}/${SAMPLENAME}.bam | grep -v ZS:i | cut -f1 | sort | uniq -c | wc -l)
MYARRAY+=("$SAMPLENAME $MYOUT")
#MYOUT=$( samtools view ${SAMPLEOUTDIR}/${SAMPLENAME}.bam | awk '$5==60' | wc -l)
#MYARRAY+=("$SAMPLENAME $MYOUT")
echo $SAMPLENAME $MYOUT
done
printf "%s\n" "${MYARRAY[@]}" > ${PROJDIR}/h2_uniq_reads_new.txt
#Count uniquely aligned reads (uniquely mapping reads have quality of 60) assuming strandedness (we try both)
module load sw/bio/samtools/0.1.18
for STRAND in RF FR
do
MYARRAY=()
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
SAMPLEOUTDIR=${PROJDIR}/h2_align_${STRAND}
MYOUT=$( samtools view ${SAMPLEOUTDIR}/${SAMPLENAME}.bam | grep -v ZS:i | cut -f1 | sort | uniq -c | wc -l)
#MYOUT=$( samtools view ${SAMPLEOUTDIR}/${SAMPLENAME}.bam | awk '$5==60' | wc -l)
MYARRAY+=("$SAMPLENAME $MYOUT")
echo $SAMPLENAME $MYOUT
done
printf "%s\n" "${MYARRAY[@]}" > ${PROJDIR}/h2_uniq_reads_${STRAND}_new.txt
done

#At the end, I found out that using strandedness we align less reads, so it was correct that the reads were unstranded!!!
#From now on, I will start again only using unstranded analysis!

#Count aligned reads
module load sw/bio/samtools/0.1.18
MYARRAY=()
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
SAMPLEOUTDIR=${PROJDIR}/h2_align
MYOUT=$( samtools view ${SAMPLEOUTDIR}/${SAMPLENAME}.bam | cut -f1 | sort | uniq -c | wc -l)
MYARRAY+=("$SAMPLENAME $MYOUT")
#MYOUT=$( samtools view ${SAMPLEOUTDIR}/${SAMPLENAME}.bam | awk '$5==60' | wc -l)
#MYARRAY+=("$SAMPLENAME $MYOUT")
echo $SAMPLENAME $MYOUT
done
printf "%s\n" "${MYARRAY[@]}" > ${PROJDIR}/h2_aligned_reads_new.txt



#Sort alignments for usage with stringtie assuming no strand
module load sw/bio/samtools/0.1.18
SAMPLEOUTDIR=${PROJDIR}/h2_align
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
samtools sort ${SAMPLEOUTDIR}/${SAMPLENAME}.bam ${SAMPLEOUTDIR}/${SAMPLENAME}.sort 
done

#Generate gtf using stringtie assuming no strand
GTFDIR=${PROJDIR}/gtf
mkdir -p $GTFDIR
SAMPLEOUTDIR=${PROJDIR}/h2_align
MYSTRINGTIE=${SOFTDIR}/stringtie-1.3.5.Linux_x86_64/stringtie
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
$MYSTRINGTIE ${SAMPLEOUTDIR}/${SAMPLENAME}.sort.bam -o ${GTFDIR}/${SAMPLENAME}.gtf
done

#Generate gtf using stringtie assuming strandedness
for STRAND in RF FR
do
GTFDIR=${PROJDIR}/gtf_${STRAND}
mkdir -p $GTFDIR
SAMPLEOUTDIR=${PROJDIR}/h2_align
MYSTRINGTIE=${SOFTDIR}/stringtie-1.3.5.Linux_x86_64/stringtie
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
$MYSTRINGTIE ${SAMPLEOUTDIR}/${SAMPLENAME}.sort.bam --{STRAND,,} -o ${GTFDIR}/${SAMPLENAME}.gtf
done
done

#Merge gtf using stringtie assuming no strand
GTFDIR=${PROJDIR}/gtf
MYSTRINGTIE=${SOFTDIR}/stringtie-1.3.5.Linux_x86_64/stringtie
$MYSTRINGTIE --merge ${GTFDIR}/LM3.gtf ${GTFDIR}/LM4.gtf ${GTFDIR}/LM6.gtf \
${GTFDIR}/LS2.gtf ${GTFDIR}/LS3.gtf ${GTFDIR}/LS6.gtf -o ${GTFDIR}/merged.gtf


#Count abundance of each gene using htseq-count assuming no strand
module load tools/htslib/latest
module load lang/python/2.7
GTFDIR=${PROJDIR}/gtf
SAMPLEOUTDIR=${PROJDIR}/h2_align
COUNTDIR=${PROJDIR}/htseq
mkdir -p $COUNTDIR
for SAMPLENAME in LM3 LM4 LM6 LS2 LS3 LS6
do
htseq-count -r pos -f bam --stranded no ${SAMPLEOUTDIR}/${SAMPLENAME}.sort.bam ${GTFDIR}/merged.gtf > ${COUNTDIR}/${SAMPLENAME}.txt 
htseq-count -r pos -f bam --stranded yes ${SAMPLEOUTDIR}/${SAMPLENAME}.sort.bam ${GTFDIR}/merged.gtf > ${COUNTDIR}/${SAMPLENAME}_yes.txt 
htseq-count -r pos -f bam --stranded reverse ${SAMPLEOUTDIR}/${SAMPLENAME}.sort.bam ${GTFDIR}/merged.gtf > ${COUNTDIR}/${SAMPLENAME}_reverse.txt 
done

#Perform DE analysis with DESeq2 using htseq-count results
module load it/lang/r/3.3
Rscript ${FUNCDIR}/functions/limodorum/00_DEseq_from_ht.r -I ${PROJDIR}/htseq -O ${PROJDIR}/DE/DE.txt
#Perform DE analysis with DESeq2 using htseq-count results using the counts separated by organism (AB and plants)
module load it/lang/r/3.3
for ORG in AB plants
do
Rscript ${FUNCDIR}/functions/limodorum/00_DEseq_from_ht.r \
-I ${PROJDIR}/htseq_${ORG} \
-O ${PROJDIR}/DE_${ORG}/${ORG}_DE_2020.txt \
-W TRUE
done

#Plot clustering of samples (all transcripts, only plants, only AB) using VST transformation
Rscript ${FUNCDIR}/functions/limodorum/00a_DEseq_plot_vst.r -I ${PROJDIR}/htseq -G ${PROJDIR}/DE/
Rscript ${FUNCDIR}/functions/limodorum/00a_DEseq_plot_vst.r -I ${PROJDIR}/htseq_AB -G ${PROJDIR}/DE_AB/
Rscript ${FUNCDIR}/functions/limodorum/00a_DEseq_plot_vst.r -I ${PROJDIR}/htseq_plants -G ${PROJDIR}/DE_plants/

