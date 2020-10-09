# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
    help="Input file of DE results with annotation [default= %default]", metavar="character"), 
  make_option(c("-D", "--defile"), type="character", default="",
    help="Input file with DE results and raw or normalized counts for the two groups [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
    help="Output file with annotations and read counts [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-I' flag.")
} else {  cat ("infile ", opt$infile, "\n")
  infile <- opt$infile  
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

add_reads<-function(infile,defile,outfile)
{
library(data.table)
library(writexl)
deres<-fread(defile,data.table=F)
cazyres<-fread(infile,data.table=F)
deres<-deres[,c("trinity_id","gene_id","LM","LS")]
cazyres<-merge(deres,cazyres,by=c("trinity_id","gene_id"),sort=F)
DEsmallcomma<-format(cazyres,decimal.mark=",")
write.table(cazyres,outfile,quote=F,row.names=F,sep="\t")
write_xlsx(DEsmallcomma,gsub(".txt",".xlsx",outfile))
}
add_reads(infile=infile,defile=defile,outfile=outfile)
