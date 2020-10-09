# Run with --help or -h flag for help.
# Written 12/31/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-G", "--GOfile"), type="character", default="",help="Formatted GO file with Phylum [default= %default]", metavar="character"), 
  make_option(c("-I", "--infile"), type="character", default="",help="De results for either fungi or plants [default= %default]", metavar="character"), 
  make_option(c("-N", "--ontology"), type="character", default="CC", help="GO ontology (can be 'BP','CC' or 'MF') [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", help="output file of GO enrichment [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$GOfile)) {
  stop("WARNING: No GOfile specified with '-G' flag.")
} else {  cat ("GOfile ", opt$GOfile, "\n")
  GOfile <- opt$GOfile  
  }

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

topgomerda<-function(x)
{
y<-x>0
return(y)
}

GO.enrich<-function(GOfile,infile,outfile,ontology,top.res=1000)
{
library(topGO)
library(data.table)
library(writexl)
myres<-fread(infile,data.table=F)
#List of all the genes obtained in the De results (specific for fungi or plants)
allgenes<-unique(unlist(strsplit(myres$Gene_id,";")))
myres$padj[is.na(myres$padj)]<-1
#Significantly DE genes
myInterestingGenes<-unique(unlist(strsplit(myres$Gene_id[myres$padj<=0.05],";")))
GO<-fread(GOfile,data.table=F)
#From the universe of GO, only extract genes to which our transcripts match
GO<-GO[GO$Gene_id%in%allgenes,]
fullset <- as.integer(allgenes %in% myInterestingGenes)
fullset<-as.numeric(fullset)
names(fullset) <- allgenes
#Create the gene2GO object
mygene2GO<-as.list(strsplit(GO$GO,";"))
mygene2GO<-setNames(mygene2GO,GO$Gene_id)
#Finally create this stupid topGO object
mytopgo<-new("topGOdata",ontology=ontology,allGenes=fullset,annot = annFUN.gene2GO, geneSel= topgomerda , gene2GO = mygene2GO)
sampleGOdata<-mytopgo
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
#Run Kolmogorov-Smirnov
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
top.res<-min(top.res,length(resultFisher@score),length(resultKS@score),length(resultKS.elim@score))
#Extract Fisher results
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
 classicKS = resultKS, elimKS = resultKS.elim,
 orderBy = "classicFisher", ranksOf = "elimKS", topNodes = top.res)
allRes<-allRes[allRes$classicFisher<=0.05,]
pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
write.table(allRes,outfile,row.names=F,quote=F,sep="\t")
allRes<-data.frame(format(allRes,decimal.mark=","))
write_xlsx(allRes,gsub(".txt",".xlsx",outfile))
}
GO.enrich(GOfile=GOfile,infile=infile,ontology=ontology,outfile=outfile)



