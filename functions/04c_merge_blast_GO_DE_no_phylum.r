# Run with --help or -h flag for help.
# Written 12/31/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-G", "--gofile"), type="character", default="",help="Formatted GO file [default= %default]", metavar="character"), 
  make_option(c("-B", "--blastfile"), type="character", default="",help="Blastx results of genes read from DE and associated with GO [default= %default]", metavar="character"), 
  make_option(c("-D", "--defile"), type="character", default="",help="DE result file of Fungi or plants [default= %default]", metavar="character"), 
  make_option(c("-O", "--out"), type="character", default="", help="output file (contains DE with GO for either fungi or plants) [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$gofile)) {
  stop("WARNING: No GO file specified with '-G' flag.")
} else {  cat ("DE file ", opt$gofile, "\n")
  gofile <- opt$gofile  
  }

if (is.null(opt$blastfile)) {
  stop("WARNING: No blast file specified with '-B' flag.")
} else {  cat ("Blast file ", opt$blastfile, "\n")
  blastfile <- opt$blastfile  
  }

if (is.null(opt$defile)) {
  stop("WARNING: No DE file specified with '-D' flag.")
} else {  cat ("DE file ", opt$defile, "\n")
  defile <- opt$defile  
  }


if (is.null(opt$out)) {
  stop("WARNING: No output file '-O' flag.")
} else {  cat ("Output file is ", opt$out, "\n")
  outfile <- opt$out  
  #setwd(wd_location)  
  }

reformat.go<-function(gofile,blastfile,defile,outfile,minfrac=0.01)
{
library(data.table)
library(writexl)
go<-fread(gofile,data.table=F)
blast<-fread(blastfile,data.table=F)
blast<-blast[,c(1,2,13)]
blast<-blast[!duplicated(blast),]
blast$Gene_id<-unlist(lapply(strsplit(unlist(lapply(strsplit(blast$V2,"\\|"),"[",4)),"\\."),"[",1))
setnames(blast,"V1","trinity_id")
blast$gene_name<-unlist(lapply(strsplit(unlist(lapply(strsplit(blast$V13,"RecName: Full="),"[",2)),";"),"[",1))
blast<-blast[,c("trinity_id","Gene_id","gene_name")]
blastgo<-merge(go,blast,by="Gene_id",all=T)
DE<-fread(defile,data.table=F)
#Start and end are the coordinates of the start and end position of the mapping of the gene reconstrcuted
#from stringtie on that reconstructed from trinity
#Sometimes we have two entries regarding that, with diffrent positions but everything else is the same.
#We treat this instances as duplicates and keep only one of them (we clean about 4000 rows out of 66000)
DE$start<-DE$end<-NULL
DE<-DE[!duplicated(DE),]
dego<-merge(DE,blastgo,by="trinity_id",all.x=T,sort=F)
DEsmall<-dego[,!names(dego)%in%c("Gene_id","GO","gene_name")]
DEsmall<-DEsmall[!duplicated(DEsmall),]
DEsmall$gene_name<-DEsmall$Gene_id<-DEsmall$GO<-""
for(aaa in 1:nrow(DEsmall))
{
	cat("Riga",aaa,"di",nrow(DEsmall),"\n")
	DEsmall$GO[aaa]<-paste(unique(sort(unlist(strsplit(dego$GO[dego$trinity_id==DEsmall$trinity_id[aaa]],";")))),sep=";",collapse=";")
	DEsmall$Gene_id[aaa]<- paste(unique(dego$Gene_id[dego$trinity_id==DEsmall$trinity_id[aaa]]),sep=";",collapse=";")
	DEsmall$gene_name[aaa]<-paste(unique(dego$gene_name[dego$trinity_id==DEsmall$trinity_id[aaa]]),sep=";",collapse=";")
}

#DEsmall<-DEsmall[!duplicated(go_small$Gene_id),]
DEsmallcomma<-format(DEsmall,decimal.mark=",")
write.table(DEsmall,outfile,quote=F,row.names=F,sep="\t")
write.table(DEsmallcomma,gsub("2020.txt","2020_comma.txt",outfile),quote=F,row.names=F,sep="\t")
write_xlsx(DEsmallcomma,gsub("2020.txt","2020.xls",outfile))
}
reformat.go(gofile=gofile,blastfile=blastfile,defile=defile,outfile=outfile)
