# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
    help="Differential expression file showing all genes (also not differentially expressed) and their functional annotation [default= %default]", metavar="character"), 
  make_option(c("-R", "--removeme"), type="character", default="description",
    help="Comma separated list of names of the columns to remove [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
    help="Output file [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")


  if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-I' flag.")
} else {  cat ("infile ", opt$infile, "\n")
  infile <- opt$infile  
  }

  if (is.null(opt$removeme)) {
  stop("WARNING: No removeme specified with '-R' flag.")
} else {  cat ("removeme ", opt$removeme, "\n")
  removeme <- opt$removeme  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

remove_by_name<-function(infile,removeme,outfile)
{
library(data.table)
library(openxlsx)
indat<-fread(infile,data.table=F)
toremove<-unlist(strsplit(removeme,","))
outdat<-indat[,!names(indat)%in%toremove]
write.table(outdat,outfile,sep="\t",quote=F,row.names=F)
newdat<-format(outdat,decimal.mark=",")
write.xlsx(newdat,gsub(".txt",".xlsx",outfile))
}
remove_by_name(infile=infile,removeme=removeme,outfile=outfile)
