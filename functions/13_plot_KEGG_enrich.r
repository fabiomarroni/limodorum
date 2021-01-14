# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-F", "--fungifile"), type="character", default="",help="KEGG enrichment file for fungi [default= %default]", metavar="character"), 
  make_option(c("-P", "--plantfile"), type="character", default="",help="KEGG enrichment file for plants [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", help="KEGG enrichment plot outfile [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$fungifile)) {
  stop("WARNING: No fungifile specified with '-F' flag.")
} else {  cat ("fungifile ", opt$fungifile, "\n")
  fungifile <- opt$fungifile  
  }

  if (is.null(opt$plantfile)) {
  stop("WARNING: No plantfile specified with '-P' flag.")
} else {  cat ("plantfile ", opt$plantfile, "\n")
  plantfile <- opt$plantfile  
  }

  
if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

plot_KEGG<-function(fungifile,plantfile,outfile,maxlength=30)
{
library(data.table)
#I will not work on fungi files because we have no significantly enriched pathways
#fungidat<-fread(fungifile,data.table=F)
plantdat<-fread(plantfile,data.table=F)
plantdat$OR[plantdat$OR=="Inf"]<-max(plantdat$OR[plantdat$OR!="Inf"])
plantdat<-plantdat[plantdat$OR>1,]
plantdat<-plantdat[plantdat$pvalue<=0.05,]
plantdat$L2FC<-log2(plantdat$OR)
#I remove data with low number of counts
plantdat<-plantdat[plantdat$Path_in_DE+plantdat$Path_in_Universe>3,]
stlen<-nchar(plantdat$KeggPath)
spacepos<-lapply(strsplit(plantdat$KeggPath," "),nchar)
len2space<-lapply(spacepos,cumsum)
whichspace<-unlist(lapply(lapply(lapply(len2space,"-",maxlength),abs),which.min))
for(aaa in 1:nrow(plantdat))
{
if(stlen[aaa]<maxlength) next
wherebreak<- whichspace[aaa]+(unlist(len2space[aaa])[whichspace[aaa]])
substr(plantdat$KeggPath[aaa],wherebreak,wherebreak)<-"\n"
}
png(outfile,width=11,height=13,units="cm",res=600,type="cairo")
par(mar=c(5,12,2,1))
barplot(plantdat$L2FC,col="dodgerblue",horiz=T,space=0.3,names.arg=plantdat$KeggPath,cex.names=0.7,las=1,main="KEGG pathway",xlab="Log 2 enrichment")
#mtext("Log 2 Enrichment",side=1,line=2,at=-3,cex=1.4)
dev.off()
browser()
}
plot_KEGG(fungifile=fungifile,plantfile=plantfile,outfile=outfile)
