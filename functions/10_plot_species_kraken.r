# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--indir"), type="character", default="",help="Input directory [default= %default]", metavar="character"), 
  make_option(c("-L", "--level"), type="character", default="G",help="Taxonomy level [default= %default]", metavar="character"),
  make_option(c("-T", "--taxid"), type="character", default="",help="List of unique taxid and their full taxonomy according to the ncbi file fullnamelineage.dmp [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", help="output plot of abundance of the selected taxonomy level [default= %default]", metavar="character")
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

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

plot_species<-function(infile,taxid,outfile, level)
{
library(data.table)
myfiles<-dir(indir,full.names=T,pattern="report.txt")
mysample<-paste("_",gsub(".kraken.report.txt","",basename(myfiles)),sep="")
finaldat<-fread(myfiles[1])
finaldat<-finaldat[finaldat$V4==level,]
finaldat<-finaldat[,c("V5","V6","V2")]
for(aaa in 2:length(myfiles))
{
secdat<-fread(myfiles[aaa])
secdat<-secdat[secdat$V4==level,]
secdat<-secdat[,c("V5","V6","V2")]
finaldat<-merge(finaldat,secdat,by=c("V5","V6"),suffixes=c(mysample[aaa-1],mysample[aaa]),all=T)
}
mytax<-scan(taxid,what="")
taxdat<-finaldat[as.character(finaldat$V5)%in%as.character(mytax),]
taxdat[is.na(taxdat)]<-0
ciccio<-grep("V2",names(taxdat))
taxdat$total<-apply(taxdat[, ..ciccio],1,sum)
taxdat<-taxdat[,c("V5","V6","total")]
taxdat<-taxdat[order(taxdat$total,decreasing=T),]
taxdat$x<-round(100*taxdat$total/sum(taxdat$total),2)
taxdat<-taxdat[1:10,]
png(outfile,width=20,height=20,units="cm",res=600,type="cairo")
par(mar=c(7,4,4,2))
barplot(taxdat$x,names.arg=taxdat$V6,las=2,ylab="Read abundance [%]")
dev.off()
}
plot_species(infile=infile,taxid=taxid,outfile=outfile,level=level)
