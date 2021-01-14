# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--indir"), type="character", default="",help="Input directory pointing to kraken output directory [default= %default]", metavar="character"), 
  make_option(c("-S", "--separate"), type="logical", default=FALSE,help="Should mychorrized and non-mychorrized plots be separated? (Only FALSE available at present) [default= %default]", metavar="character"),
  make_option(c("-L", "--level"), type="character", default="G",help="Taxonomy level [default= %default]", metavar="character"),
  make_option(c("-T", "--taxid"), type="character", default="",help="List of unique taxid and their full taxonomy according to the ncbi file fullnamelineage.dmp [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", help="output graph of fungal abundance according to kraken [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$indir)) {
  stop("WARNING: No indir specified with '-I' flag.")
} else {  cat ("indir ", opt$indir, "\n")
  indir <- opt$indir  
  }

if (is.null(opt$level)) {
  stop("WARNING: No level specified with '-L' flag.")
} else {  cat ("level ", opt$level, "\n")
  level <- opt$level  
  }

if (is.null(opt$taxid)) {
  stop("WARNING: No taxid specified with '-T' flag.")
} else {  cat ("taxid ", opt$taxid, "\n")
  taxid <- opt$taxid  
  }

if (is.null(opt$separate)) {
  stop("WARNING: No separate specified with '-S' flag.")
} else {  cat ("separate ", opt$separate, "\n")
  separate <- opt$separate  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

plot_species<-function(infile,taxid,separate,outfile, level)
{
library(data.table)
library(gplots)
allfiles<-dir(indir,full.names=T,pattern="report.txt")
# allsample<-paste("_",gsub(".kraken.report.txt","",basename(allfiles)),sep="")
if(separate==FALSE)
{
for(condition in c("LM","LS"))
{
myfiles<-allfiles[grep(condition,allfiles)]
mysample<-paste("_",gsub(".kraken.report.txt","",basename(myfiles)),sep="")
finaldat<-fread(myfiles[1])
finaldat<-finaldat[finaldat$V4==level,]
finaldat<-finaldat[,c("V5","V6","V2")]
setnames(finaldat,"V2",paste0("V2",mysample[1]))
for(aaa in 2:length(myfiles))
{
secdat<-fread(myfiles[aaa])
secdat<-secdat[secdat$V4==level,]
secdat<-secdat[,c("V5","V6","V2")]
setnames(secdat,"V2",paste0("V2",mysample[aaa]))
finaldat<-merge(finaldat,secdat,by=c("V5","V6"),all=T)
}
mytax<-scan(taxid,what="")
taxdat<-finaldat[as.character(finaldat$V5)%in%as.character(mytax),]
taxdat[is.na(taxdat)]<-0
ciccio<-grep("V2",names(taxdat))
ciccio<-grep(condition,names(taxdat))
taxdat$total<-apply(taxdat[, ..ciccio],1,sum)
taxdat$se<-apply(taxdat[, ..ciccio],1,sd)/sqrt(length(myfiles))
taxdat$ucl<-taxdat$total+1.96*taxdat$se
taxdat$lcl<-taxdat$total-1.96*taxdat$se
taxdat$use<-taxdat$total+taxdat$se
taxdat$lse<-taxdat$total-taxdat$se
taxdat<-taxdat[,c("V5","V6","total","ucl","lcl","use","lse")]
taxdat<-taxdat[order(taxdat$total,decreasing=T),]
taxdat<-taxdat[1:10,]
taxdat$x<-round(100*taxdat$total/sum(taxdat$total),2)
taxdat$ucl<-round(100*taxdat$ucl/sum(taxdat$ucl),2)
taxdat$lcl<-round(100*taxdat$lcl/sum(taxdat$lcl),2)
taxdat$use<-round(100*taxdat$use/sum(taxdat$use),2)
taxdat$lse<-round(100*taxdat$lse/sum(taxdat$lse),2)
if(condition=="LM") 
{
#These are needed to be defined only for the first plot
png(outfile,width=17,height=17,units="cm",res=600,type="cairo")
par(mar=c(8,4,4,2),mfrow=c(1,2))
}
barplot2(taxdat$x,names.arg=taxdat$V6,las=2,plot.ci=T,ci.l=taxdat$lse,ci.u=taxdat$use,ylab="Read abundance [%]",main=paste("Fungal genera in",condition))
}
dev.off()
}
}
plot_species(infile=infile,taxid=taxid,separate=separate,outfile=outfile,level=level)
