# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-F", "--fungidir"), type="character", default="",help="Folder in which GO enrichment results for fungi have been saved [default= %default]", metavar="character"), 
  make_option(c("-P", "--plantdir"), type="character", default="",help="Folder in which GO enrichment results for plants have been saved [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", help="GO enrichment plot outfile [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$fungidir)) {
  stop("WARNING: No fungidir specified with '-F' flag.")
} else {  cat ("fungidir ", opt$fungidir, "\n")
  fungidir <- opt$fungidir  
  }

  if (is.null(opt$plantdir)) {
  stop("WARNING: No plantdir specified with '-P' flag.")
} else {  cat ("plantdir ", opt$plantdir, "\n")
  plantdir <- opt$plantdir  
  }

  
if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

plot_GO<-function(fungidir,plantdir,outfile,top=10)
{
library(data.table)
fungifiles<-dir(fungidir,pattern="DE_GO_enrich_.*.txt",full.names=T)
plantfiles<-dir(plantdir,pattern="DE_GO_enrich_.*.txt",full.names=T)
for(aaa in 1:length(fungifiles))
{
fungidata<-fread(fungifiles[aaa])
plantdata<-fread(plantfiles[aaa])
#Custom cleaning of GO terms, by removing one of two closely related terms (only when they are basically the same shit) and shortening names.
#It's horrible, I know!
fungidata$Term[fungidata$Term=="positive regulation of biological proces..."]<-"positive regulation of biological process"
fungidata$Term[fungidata$Term=="process utilizing autophagic mechanism"]<-"autophagy"
fungidata<-fungidata[!duplicated(fungidata$Term),]
plantdata$Term[plantdata$Term=="organic acid transmembrane transport"]<-"carboxylic acid transmembrane transport"
plantdata$Term[plantdata$Term=="cellular potassium ion transport"]<-"potassium ion transmembrane transport"
plantdata$Term[plantdata$Term=="thiamine-containing compound metabolic p..."]<-"thiamine metabolic process"
plantdata$Term[plantdata$Term=="regulation of membrane depolarization"]<-"action potential"
plantdata$Term[plantdata$Term=="aspartic-type endopeptidase activity"]<-"aspartic-type peptidase activity"
plantdata$Term[plantdata$Term=="palmitoyl-(protein) hydrolase activity"]<-"palmitoyl hydrolase activity"
plantdata<-plantdata[!duplicated(plantdata$Term),]
fungidata<-fungidata[1:top,]
plantdata<-plantdata[1:top,]
fungidata$L2FC<-log2(fungidata$Significant/fungidata$Expected)
plantdata$L2FC<-log2(plantdata$Significant/plantdata$Expected)
if(aaa==1)
{
fungitot<-fungidata
planttot<-plantdata
}else{
fungitot<-rbind(fungitot,fungidata)
planttot<-rbind(planttot,plantdata)
}
}
mycol<-c(rep("red",top),rep("dodgerblue",top),rep("grey",top))
onto<-unlist(lapply(strsplit(basename(fungifiles),"_"),"[",4))
png(outfile,width=22,height=22,units="cm",res=600,type="cairo")
par(mar=c(4,16,2,1),mfrow=c(1,2),xpd=NA)
barplot(fungitot$L2FC,col=mycol,horiz=T,names.arg=fungitot$Term,cex.names=0.9,las=1,main="Fungi")
barplot(planttot$L2FC,col=mycol,horiz=T,names.arg=planttot$Term,cex.names=0.9,las=1,main="Plants")
#legend("bottomright",legend=rev(onto),fill=rev(unique(mycol)),bty="n")
legend(x=1.6,y=3,legend=rev(onto),fill=rev(unique(mycol)),bty="n")
mtext("Log 2 Enrichment",side=1,line=2,at=-6,cex=1.2)
dev.off()
}
plot_GO(fungidir=fungidir,plantdir=plantdir,outfile=outfile)
