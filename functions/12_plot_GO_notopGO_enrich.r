# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-F", "--fungidir"), type="character", default="",help="Folder in which GO enrichment results for fungi have been saved [default= %default]", metavar="character"), 
  make_option(c("-P", "--plantdir"), type="character", default="",help="Folder in which GO enrichment results for plants have been saved [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", help="GO enrichment (based on simple Fisher test, NOT in topGO implementation) plot outfile [default= %default]", metavar="character"),
  make_option(c("-D", "--outdir"), type="character", default="", help="GO enrichment table directory [default= %default]", metavar="character")
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

  if (is.null(opt$outdir)) {
  stop("WARNING: No outdir specified with '-D' flag.")
} else {  cat ("outdir ", opt$outdir, "\n")
  outdir <- opt$outdir  
  }

 
if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

plot_GO<-function(fungidir,plantdir,outfile,outdir,top=10)
{
library(data.table)
fungifiles<-dir(fungidir,pattern="DE_GO_notopGO_enrich_.*.txt",full.names=T)
plantfiles<-dir(plantdir,pattern="DE_GO_notopGO_enrich_.*.txt",full.names=T)
for(aaa in 1:length(fungifiles))
{
fungidata<-fread(fungifiles[aaa])
plantdata<-fread(plantfiles[aaa])

fungidata<-fungidata[1:top,]
plantdata<-plantdata[1:top,]
#######Custom changes to beautify GO names
fungidata$GOterm[fungidata$GOterm=="sporulation resulting in formation of a cellular spore"]<-"sporulation (cellular spore)"
plantdata$GOterm[plantdata$GOterm=="oxidoreductase activity, acting on NAD(P)H, oxygen as acceptor"]<-"oxidoreductase activity, acting on NAD(P)H"
plantdata$GOterm[plantdata$GOterm=="positive regulation of cell communication by electrical coupling"]<-"positive regulation of cell communication..."
plantdata$GOterm[plantdata$GOterm=="proteolysis involved in cellular protein catabolic process"]<-"proteolysis involved in (...) catabolic process"
plantdata$GOterm[plantdata$GOterm=="programmed cell death involved in cell development"]<-"programmed cell death..."
plantdata$GOterm[plantdata$GOterm=="hydrolase activity, hydrolyzing O-glycosyl compounds"]<- "hydrolase activity (O-glycosyl compounds)"

fungidata$OR[fungidata$OR=="Inf"]<-max(fungidata$OR[fungidata$OR!="Inf"])
plantdata$OR[plantdata$OR=="Inf"]<-max(plantdata$OR[plantdata$OR!="Inf"])
fungidata<-fungidata[order(fungidata$OR,decreasing=F),]
plantdata<-plantdata[order(plantdata$OR,decreasing=F),]
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
onto<-unlist(lapply(strsplit(basename(fungifiles),"_"),"[",5))
png(outfile,width=20,height=20,units="cm",res=600,type="cairo")
par(mar=c(4,15,2,1),mfrow=c(1,2),xpd=NA)
barplot(fungitot$OR,col=mycol,horiz=T,names.arg=fungitot$GOterm,cex.names=0.9,las=1,main="Fungi")
barplot(planttot$OR,col=mycol,horiz=T,names.arg=planttot$GOterm,cex.names=0.9,las=1,main="Plants")
#legend("bottomright",legend=rev(onto),fill=rev(unique(mycol)),bty="n")
legend(x=1.9,y=3,legend=rev(onto),fill=rev(unique(mycol)),bty="n")
mtext("Odds Ratio",side=1,line=2,at=-45,cex=1.2)
dev.off()
write.table(fungitot,paste(outdir,"AB_enrichnotopGO.txt",sep="/"),quote=F,row.names=F,sep="\t")
write.table(planttot,paste(outdir,"plants_enrichnotopGO.txt",sep="/"),quote=F,row.names=F,sep="\t")
}
plot_GO(fungidir=fungidir,plantdir=plantdir,outdir=outdir,outfile=outfile)
