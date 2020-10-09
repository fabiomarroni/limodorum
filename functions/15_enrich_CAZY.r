# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-P", "--plantfile"), type="character", default="",
    help="Differential expression file showing all genes (also not differentially expressed) and their functional annotation in plants[default= %default]", metavar="character"), 
  make_option(c("-F", "--fungifile"), type="character", default="",
    help="Differential expression file showing all genes (also not differentially expressed) and their functional annotation in fungi [default= %default]", metavar="character"), 
  make_option(c("-M", "--mintools"), type="numeric", default=2,
    help="Minimum number of tools assigning a CAZY for result to be considered as positive [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
    help="CAZy enrichment output graph [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")


  if (is.null(opt$plantfile)) {
  stop("WARNING: No plantfile specified with '-P' flag.")
} else {  cat ("plantfile ", opt$plantfile, "\n")
  plantfile <- opt$plantfile  
  }

  if (is.null(opt$fungifile)) {
  stop("WARNING: No fungifile specified with '-F' flag.")
} else {  cat ("fungifile ", opt$fungifile, "\n")
  fungifile <- opt$fungifile  
  }

  if (is.null(opt$mintools)) {
  stop("WARNING: No mintools specified with '-M' flag.")
} else {  cat ("mintools ", opt$mintools, "\n")
  mintools <- opt$mintools  
  }

  
if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

enrich_cazy<-function(plantfile,fungifile,mintools,outfile)
{
library(data.table)
png(outfile,width=12,height=10,units="cm",res=600,type="cairo")
par(mar=c(4,5,2,1),mfrow=c(1,2))
for(defile in c(fungifile,plantfile))
{
mymain<-ifelse(defile==fungifile,"Fungi","Plants")
deres<-fread(defile,data.table=F)
deres$any<-apply(deres[,c("HMMER","Hotpep","DIAMOND")],1,paste,collapse=";")
deres<-deres[,c("trinity_id","padj","any","#ofTools")]
sigpos<-na.omit(deres[deres$padj<=0.05&deres$"#ofTools">=mintools,])
signeg<-na.omit(deres[deres$padj<=0.05&deres$"#ofTools"<mintools,])
nosigpos<-deres[(is.na(deres$padj)|deres$padj>0.05)&deres$"#ofTools">=mintools,]
nosigneg<-deres[(is.na(deres$padj)|deres$padj>0.05)&deres$"#ofTools"<mintools,]

cazyClasses<-rev(c("GH","GT","PL","CE","AA"))
spcount<-rep(0,length(cazyClasses))
sncount<-rep(0,length(cazyClasses))
npcount<-rep(0,length(cazyClasses))
nncount<-rep(0,length(cazyClasses))
fisherp<-OR<-rep(0,length(cazyClasses))

names(spcount)<-names(sncount)<-names(npcount)<-names(nncount)<-names(OR)<-names(fisherp)<-cazyClasses
for(aaa in 1:length(cazyClasses))
{
spcount[aaa]<-length(grep(cazyClasses[aaa],sigpos$any))
#For each Cazy class the negatives are those that are associated to a DIFFERENT cazy class (spcount that do not grep with the class)
#plus those that are negative for ALL classses (nrow of negatives)
sncount[aaa]<-length(grep(cazyClasses[aaa],sigpos$any,invert=TRUE))+nrow(signeg)
npcount[aaa]<-length(grep(cazyClasses[aaa],nosigpos$any))
#For each Cazy class the negatives are those that are associated to a DIFFERENT cazy class (spcount that do not grep with the class)
#plus those that are negative for ALL classses (nrow of negatives)
nncount[aaa]<-length(grep(cazyClasses[aaa],nosigpos$any,invert=TRUE))+nrow(nosigneg)
fres<-fisher.test(matrix(c(spcount[aaa],sncount[aaa],npcount[aaa],nncount[aaa]),nrow=2))
fisherp[aaa]<-fres$p.value
OR[aaa]<-fres$estimate
}
myast<-rep("",length(cazyClasses))
myast[fisherp<=0.05]<-"*"
mynames<-ifelse(defile==fungifile,cazyClasses,rep("",length(cazyClasses)))
onlydata<-barplot(OR,xlim=c(0,4),col="dodgerblue",horiz=T,cex.names=0.9,las=1,main=mymain,plot=F)
barplot(OR,xlim=c(0,4.5),col="dodgerblue",horiz=T,names.arg=cazyClasses,cex.names=0.9,las=1,main=mymain)
text(x=4.2,y=onlydata,labels=myast,cex=2.5)
}
mtext("Odds Ratio",side=1,line=3,at=-3.5,cex=1.4)
dev.off()
}
enrich_cazy(plantfile=plantfile,fungifile=fungifile,mintools=mintools,outfile=outfile)
