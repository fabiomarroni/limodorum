# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",help="DE file with KEGG annotation [default= %default]", metavar="character"), 
  make_option(c("-R", "--regulation"), type="character", default="",help="Direction of regulation. Possible values are 'up', 'down', 'all', and 'no' [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", help="Output file of KEGG enrichment analysis [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-I' flag.")
} else {  cat ("Infile ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$regulation)) {
  stop("WARNING: No regulation specified with '-R' flag.")
} else {  cat ("Regulation is", opt$regulation, "\n")
  regulation <- opt$regulation  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

kegg.enrich<-function(infile,regulation,outfile)
{
library(data.table)
library(writexl)
indata<-fread(infile,data.table=F)
myKP<-unique(unlist(strsplit(indata$Kpathway,";")))
myres<-data.frame(KeggPath=myKP,Path_in_DE=NA,No_path_in_DE=NA,Path_in_Universe=NA,No_path_in_Universe=NA,OR=NA,pvalue=NA)
indata$padj[is.na(indata$padj)]<-1
#Define the set of interest (DE) and the Universe based on the kind of regulation we want to study
if(regulation=="all")
{
	DE<-indata[indata$padj<=0.05,]
	Universe<-indata[indata$padj>0.05,]
}
if(regulation=="up")
{
	DE<-indata[indata$padj<=0.05&indata$log2FoldChange>0,]
	Universe<-indata[indata$padj>0.05|indata$log2FoldChange<=0,]
}
if(regulation=="down")
{
	DE<-indata[indata$padj<=0.05&indata$log2FoldChange<0,]
	Universe<-indata[indata$padj>0.05|indata$log2FoldChange>=0,]
}
#In this case we just want to see enrichment in not significantly DE genes
if(regulation=="no")
{
	DE<-indata[indata$padj>0.05,]
	Universe<-indata[indata$padj<=0.05,]
}
for(aaa in 1:nrow(myres))
{
myres$Path_in_DE[aaa]<-length(grep(myres$KeggPath[aaa],DE$Kpathway))
myres$No_path_in_DE[aaa]<-length(grep(myres$KeggPath[aaa],DE$Kpathway,invert=TRUE))
myres$Path_in_Universe[aaa]<-length(grep(myres$KeggPath[aaa],Universe$Kpathway))
myres$No_path_in_Universe[aaa]<-length(grep(myres$KeggPath[aaa],Universe$Kpathway,invert=TRUE))
fr<-fisher.test(matrix(c(myres$Path_in_DE[aaa],myres$No_path_in_DE[aaa],myres$Path_in_Universe[aaa],myres$No_path_in_Universe[aaa]),
		nrow=2,byrow=T))
myres$pvalue[aaa]<-fr$p.value
myres$OR[aaa]<-fr$estimate
}
myres<-myres[order(myres$pvalue),]
write.table(myres,outfile,sep="\t",row.names=F,quote=F)
myres<-data.frame(format(myres,decimal.mark=","))
write_xlsx(myres,gsub(".txt",".xlsx",outfile))
}
kegg.enrich(infile=infile,regulation=regulation,outfile=outfile)
