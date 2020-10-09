

# Run with --help or -h flag for help.
# Written 12/31/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",help="File with DE results and annotaiton for plants or fungi [default= %default]", metavar="character"), 
  make_option(c("-N", "--ontology"), type="character", default="CC", help="Go ontology (can be 'BP','CC' or 'MF') [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", help="GO enrichment results computed with Fisher's test NOT in topGO implementation [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-I' flag.")
} else {  cat ("Infile ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$ontology)) {
  stop("WARNING: No ontology specified with '-N' flag.")
} else {  cat ("ontology ", opt$ontology, "\n")
  ontology <- opt$ontology  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }


GO.enrich<-function(infile,outfile,ontology,top.res=1000)
{
library(data.table)
library(writexl)
library(GO.db)
indata<-fread(infile,data.table=F)
allGO<-unique(unlist(strsplit(indata$GO,";")))
myonto<-GOTERM[Ontology(GOTERM)==ontology]
#Only select GO terms in the current ontology
myGO<-GOID(myonto)[names(GOID(myonto))%in%allGO]
myres<-data.frame(GO=myGO,GO_in_DE=NA,No_GO_in_DE=NA,GO_in_Universe=NA,No_GO_in_Universe=NA,OR=NA,pvalue=NA,GOterm=NA,stringsAsFactors=F)
indata$padj[is.na(indata$padj)]<-1
DE<-indata[indata$padj<=0.05,]
Universe<-indata[indata$padj>0.05,]
for(aaa in 1:nrow(myres))
{
cat("Row",aaa,"out of",nrow(myres),"\n")
myres$GO_in_DE[aaa]<-length(grep(myres$GO[aaa],DE$GO))
myres$No_GO_in_DE[aaa]<-length(grep(myres$GO[aaa],DE$GO,invert=TRUE))
myres$GO_in_Universe[aaa]<-length(grep(myres$GO[aaa],Universe$GO))
myres$No_GO_in_Universe[aaa]<-length(grep(myres$GO[aaa],Universe$GO,invert=TRUE))
fr<-fisher.test(matrix(c(myres$GO_in_DE[aaa],myres$No_GO_in_DE[aaa],myres$GO_in_Universe[aaa],myres$No_GO_in_Universe[aaa]),
		nrow=2,byrow=T))
myres$pvalue[aaa]<-fr$p.value
myres$OR[aaa]<-fr$estimate
myres$GOterm[aaa]<-Term(myonto)[names(Term(myonto))==myres$GO[aaa]]
}

myres<-myres[order(myres$pvalue),]
write.table(myres,outfile,row.names=F,quote=F,sep="\t")
}
GO.enrich(infile=infile,ontology=ontology,outfile=outfile)



