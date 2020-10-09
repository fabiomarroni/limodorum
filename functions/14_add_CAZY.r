# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--indir"), type="character", default="",
    help="Directory in which raw CAZY results files are stored [default= %default]", metavar="character"), 
  make_option(c("-D", "--defile"), type="character", default="",
    help="Differential expression file showing all genes (also not differentially expressed) and their functional annotation [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
    help="Output file, equal to 'defile' with the addition of CAZY results [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$indir)) {
  stop("WARNING: No indir specified with '-I' flag.")
} else {  cat ("indir ", opt$indir, "\n")
  indir <- opt$indir  
  }

  if (is.null(opt$defile)) {
  stop("WARNING: No defile specified with '-D' flag.")
} else {  cat ("defile ", opt$defile, "\n")
  defile <- opt$defile  
  }

  
if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

add_cazy<-function(indir,defile,outfile)
{
library(data.table)
library(writexl)
deres<-fread(defile,data.table=F)
tomerge<-dir(indir,pattern="out.txt",full.names=T)
fullres<-fread(tomerge[1],data.table=F)
#Loop over the files and merge them, after checking that they differ (they can be the same if I wrote twice the same results to different files)
for(aaa in 2:length(tomerge))
{
tres<-fread(tomerge[aaa],data.table=F)
if(sum(tres$"Gene ID"%in%fullres$"Gene ID")>nrow(tres)/2) stop("Error! Too many duplicate genes observed in file", tomerge[aaa])
fullres<-rbind(fullres,tres)
}
#Since dbCAN returns transcripts broken down by (supposed) functional domains, I merge them in a single transcript.
#The different assignments are then reported in the same line separated by a semi-colon
tlist<-lapply(strsplit(fullres$"Gene ID","_"),"[",1:5)
fullres$gene_id<-unlist(lapply(tlist,paste,sep="_",collapse="_"))
#fullres[fullres=="N"]<-""
newres<-data.frame(gene_id=unique(fullres$gene_id),stringsAsFactors=F)
newres$"#ofTools"<-newres$Signalp<-newres$DIAMOND<-newres$Hotpep<-newres$HMMER<-NA
for(bbb in 1:nrow(newres))
{
newres$HMMER[bbb]<-paste(fullres$HMMER[fullres$gene_id==newres$gene_id[bbb]],sep=";",collapse=";")
newres$Hotpep[bbb]<-paste(fullres$Hotpep[fullres$gene_id==newres$gene_id[bbb]],sep=";",collapse=";")
newres$DIAMOND[bbb]<-paste(fullres$DIAMOND[fullres$gene_id==newres$gene_id[bbb]],sep=";",collapse=";")
newres$Signalp[bbb]<-paste(fullres$Signalp[fullres$gene_id==newres$gene_id[bbb]],sep=";",collapse=";")
newres$"#ofTools"[bbb]<-paste(fullres$"#ofTools"[fullres$gene_id==newres$gene_id[bbb]],sep=";",collapse=";")
}
finderes<-merge(deres,newres,by.x="trinity_id",by.y="gene_id",all.x=T)
finderes$HMMER[is.na(finderes$HMMER)]<-""
finderes$Hotpep[is.na(finderes$Hotpep)]<-""
finderes$DIAMOND[is.na(finderes$DIAMOND)]<-""
finderes$Signalp[is.na(finderes$Signalp)]<-""
finderes$"#ofTools"[is.na(finderes$"#ofTools")]<-""
finderes<-finderes[order(finderes$padj),]
DEsmallcomma<-format(finderes,decimal.mark=",")
write.table(finderes,outfile,quote=F,row.names=F,sep="\t")
write_xlsx(DEsmallcomma,gsub(".txt",".xlsx",outfile))
}
add_cazy(indir=indir,defile=defile,outfile=outfile)
