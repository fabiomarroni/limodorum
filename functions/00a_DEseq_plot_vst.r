# Run with --help flag for help.
# Modified 12/30/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--input_dir"), type="character", default="",
              help="Input directory with htseq-count output files[default= %default]", metavar="character"),
  make_option(c("-G", "--graphdir"), type="character", default="", 
              help="output graph directory [default= %default]", metavar="character")
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

if (is.null(opt$graphdir)) {
  stop("WARNING: No graphdir directory specified with '-I' flag.")
} else {  cat ("graphdir is ", opt$graphdir, "\n")
  graphdir <- opt$graphdir  
  }

runDEseq<-function(starcount, graphdir) 
{
    library("pheatmap")
    library("ggplot2")
    library("ggrepel")
	library(DESeq2)
    library(data.table)
    sampleFiles <- list.files(starcount)
	#Ignore the total statistics file
	sampleFiles <- sampleFiles[!sampleFiles%in%"total.txt"]    
    sampleCondition <- substr(sampleFiles,1,2)
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
	#Plot some clustering
	vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=FALSE)
	select <- order(rowMeans(counts(ddsHTSeq, normalized=TRUE)), decreasing=TRUE)[1:100]
	pheatmap(assay(vsd)[select[1:100],], cluster_rows=FALSE, fontsize=4, cellwidth=6, cellheight=4, filename=paste(graphdir,"100-most-expressed-genes_clust.pdf",sep=""))
	#dev.off()
	mydat<-data.frame(assay(vsd),check.names=F)
	myvar<-apply(mydat,1,var)
	mydat<-mydat[order(myvar,decreasing=T),]
	mydat<-mydat[1:500,]
	mypc<-prcomp(t(mydat))
	percentVar <- round(100*mypc$sdev^2/sum(mypc$sdev^2),0)
	myxlab<-paste("PC1: ",percentVar[1],"% variance",sep="")
	myylab<-paste("PC2: ",percentVar[2],"% variance",sep="")
    pdf(paste(graphdir,"VST_PCA.pdf",sep=""))
	pino<-plotPCA(vsd, intgroup=c("condition"),returnData=T)
	pino$name<-gsub(".txt","",pino$name)
	myplot<-ggplot(pino, aes(x=PC1, y=PC2, color=group),size=3) +geom_point(size=3) +geom_text_repel(aes(label=name, color=group),size=4) + coord_fixed() + 	theme(legend.text=element_text(size=10)) +xlab(myxlab) + ylab(myylab)
	print(myplot)	
	dev.off()
}
runDEseq(starcount=starcount,graphdir=graphdir)
