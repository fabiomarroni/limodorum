# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-G", "--GOdir"), type="character", default="",help="Folder in which GO enrichment results have been saved [default= %default]", metavar="character"), 
  make_option(c("-K", "--keggfile"), type="character", default="",help="File in which KEGG enrichment results for plants have been saved [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", help="GO and KEGG enrichment plot outfile [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$GOdir)) {
  stop("WARNING: No GOdir specified with '-F' flag.")
} else {  cat ("GOdir ", opt$GOdir, "\n")
  GOdir <- opt$GOdir  
  }

  if (is.null(opt$keggfile)) {
  stop("WARNING: No keggfile specified with '-P' flag.")
} else {  cat ("keggfile ", opt$keggfile, "\n")
  keggfile <- opt$keggfile  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

plot_GO<-function(GOdir,keggfile,outfile,top=10,maxlength=25)
{
library(data.table)


######################################
#PLOT GO PANEL
######################################
plantfiles<-dir(GOdir,pattern="DE_GO_notopGO_enrich_.*.txt",full.names=T)
for(aaa in 1:length(plantfiles))
{
plantdata<-fread(plantfiles[aaa])

plantdata<-plantdata[1:top,]
#######Custom changes to beautify GO names
plantdata$GOterm[plantdata$GOterm=="oxidoreductase activity, acting on NAD(P)H, oxygen as acceptor"]<-"oxidoreductase activity, acting on NAD(P)H"
plantdata$GOterm[plantdata$GOterm=="positive regulation of cell communication by electrical coupling"]<-"positive regulation of cell communication..."
plantdata$GOterm[plantdata$GOterm=="proteolysis involved in cellular protein catabolic process"]<-"proteolysis involved in (...) catabolic process"
plantdata$GOterm[plantdata$GOterm=="programmed cell death involved in cell development"]<-"programmed cell death..."
plantdata$GOterm[plantdata$GOterm=="hydrolase activity, hydrolyzing O-glycosyl compounds"]<- "hydrolase activity (O-glycosyl compounds)"

plantdata$OR[plantdata$OR=="Inf"]<-max(plantdata$OR[plantdata$OR!="Inf"])
plantdata$OR[plantdata$OR==0]<-min(plantdata$OR[plantdata$OR!=0])
plantdata<-plantdata[order(plantdata$OR,decreasing=F),]
if(aaa==1)
{
planttot<-plantdata
}else{
planttot<-rbind(planttot,plantdata)
}
}
planttot$log2enrich<-log(planttot$OR,2)
planttot<-planttot[planttot$log2enrich>0,]
mycol<-c(rep("red",top),rep("dodgerblue",top),rep("grey",top))
onto<-unlist(lapply(strsplit(basename(plantfiles),"_"),"[",5))
#browser()
png(outfile,width=20,height=15,units="cm",res=600,type="cairo")
par(mar=c(4,15,2,0),mfrow=c(1,2),xpd=NA)
barplot(planttot$log2enrich,col=mycol,horiz=T,names.arg=planttot$GOterm,cex.names=0.8,las=1,main="GO")
#legend("bottomright",legend=rev(onto),fill=rev(unique(mycol)),bty="n")
legend(x=2.1,y=4.2,legend=rev(onto),fill=rev(unique(mycol)),bty="n")
#mtext("Odds Ratio",side=1,line=2,at=-45,cex=1.2)



######################################
#PLOT KEGG PANEL
######################################
plantdat<-fread(keggfile,data.table=F)
plantdat<-plantdat[plantdat$OR>1,]
plantdat<-plantdat[plantdat$pvalue<=0.05,]
plantdat$L2FC<-log2(plantdat$OR)
stlen<-nchar(plantdat$KeggPath)
spacepos<-lapply(strsplit(plantdat$KeggPath," "),nchar)
len2space<-lapply(spacepos,cumsum)
whichspace<-unlist(lapply(lapply(lapply(len2space,"-",maxlength),abs),which.min))
for(aaa in 1:nrow(plantdat))
{
if(stlen[aaa]<maxlength) next
wherebreak<- whichspace[aaa]+(unlist(len2space[aaa])[whichspace[aaa]])
substr(plantdat$KeggPath[aaa],wherebreak,wherebreak)<-"\n"
}
par(mar=c(4,12,2,1),xpd=NA)
barplot(plantdat$L2FC,col="dodgerblue",horiz=T,names.arg=plantdat$KeggPath,cex.names=0.8,las=1,main="KEGG")
mtext("Log 2 enrichment",side=1,line=2,at=-5,cex=1.4)

dev.off()
}
plot_GO(GOdir=GOdir,keggfile=keggfile,outfile=outfile)
