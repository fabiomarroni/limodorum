# Run with --help or -h flag for help.
# Written 12/31/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-D", "--DE_species"), type="character", default="",help="Result of DEseq to which we attached the species and taxid using Fabio's script [default= %default]", metavar="character"), 
  make_option(c("-T", "--taxid"), type="character", default="",help="List of unique taxid and their full taxonomy according to the ncbi file fullnamelineage.dmp [default= %default]", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="", help="output file (contains DE results and assignment to Phyla) [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$DE_species)) {
  stop("WARNING: No (modified) DE file '-D' flag.")
} else {  cat ("DE file ", opt$DE_species, "\n")
  DE_species <- opt$DE_species  
  }

if (is.null(opt$taxid)) {
  stop("WARNING: No file with taxonomy of transcripts in DE file '-T' flag.")
} else {  cat ("Taxonomy file ", opt$taxid, "\n")
  taxid <- opt$taxid  
  }

if (is.null(opt$out)) {
  stop("WARNING: No output file '-O' flag.")
} else {  cat ("Output file is ", opt$out, "\n")
  outfile <- opt$out  
  #setwd(wd_location)  
  }

sep_plants_fungi<-function(DE_species, taxid, outfile) 
{
    library(data.table)
	library(writexl)
	DEspecies<-fread(DE_species,data.table=F,fill=T)
	DEspecies$Taxid<-unlist(lapply(strsplit(DEspecies$description,"taxid_"),"[",2))	
	DEspecies$Taxid<-gsub("\\)","",DEspecies$Taxid)
	taxdat<-fread(taxid,data.table=F,fill=TRUE)
	DEnew<-DEspecies[DEspecies$Taxid%in%unlist(taxdat),]	
	DEnew$Taxid<-NULL	
	write.table(DEnew,outfile,sep="\t",quote=F,row.names=F)
	# newresults<-format(DEnew,decimal.mark=",")
	# write_xlsx(newresults,gsub(".txt",".xlsx",outfile))
}
sep_plants_fungi(DE_species, taxid=taxid, outfile)
