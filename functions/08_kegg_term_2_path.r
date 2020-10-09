# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",help="DE file with KEGG annotation [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", help="output file (contains DE with GO, KEGG and KEGG pathway annotation) [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-I' flag.")
} else {  cat ("Infile ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

assign.kegg.pathway<-function(infile,outfile)
{
library(data.table)
library(KEGGREST)
library(writexl)
indata<-fread(infile,data.table=F)
indata$gene_name<-""
indata$KO<-""
indata$Kpathway<-""
for(aaa in 1:nrow(indata))
{
    cat("Record", aaa, "out of",nrow(indata),"\n")	
	#if(aaa==416) browser()
    if(indata$KEGG[aaa]=="") next
  	#We sometimes have two kegg terms together, separated by ";". 
	#We have to split them and search for them separately. We put them in the vector
    #fullkegg and loop over it
    fullkegg<-unique(unlist(strsplit(indata$KEGG[aaa],";")))
	fullkegg<-fullkegg[fullkegg!=""]    
	if(length(fullkegg)==0) next
	myname<-myorthology<-mypathway<-rep("",length(fullkegg))
	for(bbb in 1:length(fullkegg))
    {
	KI<-tryCatch(keggGet(fullkegg[bbb]),error = function(e) {NULL})
    if(!is.null(KI[[1]]$NAME)) myname[bbb]<-KI[[1]]$NAME
    if(!is.null(KI[[1]]$ORTHOLOGY)) myorthology[bbb]<-KI[[1]]$ORTHOLOGY
    if(!is.null(KI[[1]]$PATHWAY)) mypathway[bbb]<-paste(KI[[1]]$PATHWAY,collapse=";")
    # if(!is.null(KI[[1]]$NAME)) indata$gene_name[aaa]<-KI[[1]]$NAME
    # if(!is.null(KI[[1]]$ORTHOLOGY)) indata$KO[aaa]<-KI[[1]]$ORTHOLOGY
    # if(!is.null(KI[[1]]$PATHWAY)) indata$Kpathway[aaa]<-KI[[1]]$PATHWAY
    }
    indata$gene_name[aaa]<-paste(unique(myname),collapse=";")
    indata$KO[aaa]<-paste(unique(myorthology),collapse=";")
    indata$Kpathway[aaa]<-paste(unique(unlist(strsplit(mypathway,";"))),collapse=";")
}
write.table(indata,outfile,sep="\t",row.names=F,quote=F)
allRes<-data.frame(format(indata,decimal.mark=","))
write_xlsx(allRes,gsub(".txt",".xlsx",outfile))

}
assign.kegg.pathway(infile=infile,outfile=outfile)
