# Run with --help flag for help.
# Modified 12/30/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--input_dir"), type="character", default="",
              help="Input directory containing htseq-count output files [default= %default]", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="", 
              help="DE output file name [default= %default]", metavar="character"),
  make_option(c("-W", "--write.counts"), type="logical", default=TRUE, 
              help="should count tables be written? [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$input)) {
  stop("WARNING: No input directory specified with '-I' flag.")
} else {  cat ("Working directory is ", opt$input, "\n")
  starcount <- opt$input_dir  
  #setwd(wd_location)  
  }

if (is.null(opt$out)) {
  stop("WARNING: No input directory specified with '-I' flag.")
} else {  cat ("Output file is ", opt$out, "\n")
  outfile <- opt$out  
  #setwd(wd_location)  
  }

if (is.null(opt$write.counts)) {
  stop("WARNING: No write.counts specified with '-W' flag.")
} else {  cat ("write.counts is ", opt$write.counts, "\n")
  write.counts <- opt$write.counts  
  }

runDEseq<-function(starcount, outfile,write.counts) 
{
    library(DESeq2)
    library(data.table)
    library(openxlsx)
	sampleFiles <- list.files(starcount)
	#Ignore the total statistics file
	sampleFiles <- sampleFiles[!sampleFiles%in%"total.txt"]    
	sampleCondition <- rev(substr(sampleFiles,1,2))
	sampleTable <- data.frame(sampleName = sampleFiles,
                              fileName = sampleFiles,
                              condition = sampleCondition)
    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                        directory = starcount,
                                        design= ~ condition)
	#Remove lowly expressed genes
	keep <- rowSums(counts(ddsHTSeq)) >= 10    
	ddsHTSeq<-ddsHTSeq[keep,]
	ddsHTSeq<-DESeq(ddsHTSeq)
	res <- results(ddsHTSeq)
	res<-res[order(res$padj),]
	cat("Writing results\n")    
	write.table(res,outfile,sep="\t",quote=F)
	if(write.counts)
	{
	myraw<-counts(ddsHTSeq,normalized=F)
	myrawLM<-myraw[,grep("LM",colnames(myraw))]
	myrawLS<-myraw[,grep("LS",colnames(myraw))]
	nraw<-cbind(apply(myrawLM,1,mean),apply(myrawLS,1,mean))
	nraw<-round(nraw,0)
	colnames(nraw)<-c("LM","LS")
	raw.res<-merge(res,nraw,by="row.names",sort=F,all=T)
	write.table(raw.res,gsub("_2020.txt","_raw_2020.txt",outfile),sep="\t",quote=F)
	#We still call them raw, but we use normalized!
	myraw<-counts(ddsHTSeq,normalized=T)
	myrawLM<-myraw[,grep("LM",colnames(myraw))]
	myrawLS<-myraw[,grep("LS",colnames(myraw))]
	nraw<-cbind(apply(myrawLM,1,mean),apply(myrawLS,1,mean))
	nraw<-round(nraw,0)
	colnames(nraw)<-c("LM","LS")
	raw.res<-merge(res,nraw,by="row.names",sort=F,all=T)
	write.table(raw.res,gsub("_2020.txt","_norm_2020.txt",outfile),sep="\t",quote=F)
	}
}
runDEseq(starcount=starcount,outfile=outfile,write.counts=write.counts)
