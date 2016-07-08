from __future__ import with_statement 

# ==============================================================================
# 						plotMyGBK
#
# Author: Sandro Valenzuela (sandrolvalenzuead@gmail.com) 
#
# Please type "python plotMyGBK.py -h" for usage help
#
# ==============================================================================

__author__ = 'Sandro Valenzuela (sandrolvalenzuead@gmail.com)'
__version__ = '1.0'
__date__ = '15 February 2016'

import sys, os, re, multiprocessing
import subprocess, csv
from optparse import OptionParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

def printPlotStep():
	plotstep=open("plotstep.R", 'w')
	plotstep.write("""rm(list=ls());
library(Rsamtools)
library(OmicCircos)
library(data.table)
args<-commandArgs()
gl<-as.numeric(args[6])
makecog<-c(args[7])
measure<-data.frame(seg.name=seq(1:8), seg.start=seq(1:8), seg.end=seq(1:8)+1, seg.value=1)
measure2<-data.frame(seg.name=1, seg.start=1, seg.end=signif(2*gl/12/1000), seg.value=round(2*gl/12/1000,0))
measure4<-data.frame(seg.name=1, seg.start=1, seg.end=signif(4*gl/12/1000), seg.value=round(4*gl/12/1000,0))
measure8<-data.frame(seg.name=1, seg.start=1, seg.end=signif(8*gl/12/1000), seg.value=round(8*gl/12/1000,0))
measure10<-data.frame(seg.name=1, seg.start=1, seg.end=signif(10*gl/12/1000), seg.value=round(10*gl/12/1000,0))
measure["V5"]<-measure2["V5"]<-measure4["V5"]<-measure8["V5"]<-measure10["V5"]<-1
data<-read.table("contigplot.dat",header=F)
fdata<-read.table("forwardplot.dat",header=F)
rdata<-read.table("reverseplot.dat",header=F)
rnadata<-read.table("rna.dat",header=F)
gcdata<-read.table("gcskewplot.dat",header=F)
if(nrow(data)>1){
  contig<-gcdata["V1"][1,]
  sum<-0
  positive<-data.frame(name=1, value=1)
  negative<-data.frame(name=1, value=1)
  pcont=1
  ncont=1
  for(i in seq(1:nrow(gcdata))){
    if(contig==gcdata["V1"][i,]){
      sum<-sum+gcdata["V5"][i,]
    }else{
      if(sum>0){
        positive[pcont,1]<-as.character(contig)
        positive[pcont,2]<-sum
        pcont<-pcont+1
      }else{
        negative[ncont,1]<-as.character(contig)
        negative[ncont,2]<-sum
        ncont<-ncont+1
      }
      contig<-gcdata["V1"][i,]
      sum<-0
      sum<-sum+gcdata["V5"][i,]
    }
  }
  if(sum>0){
    positive[pcont,1]<-as.character(contig)
    positive[pcont,2]<-sum
  }else{
    negative[ncont,1]<-as.character(contig)
    negative[ncont,2]<-sum
  }
  positive<-positive[order(positive$value),]
  negative<-negative[order(negative$value),]
  if(nrow(positive>1)){
    if(nrow(negative)>1){
      newcontig<-rbind(positive,negative)
    }else{
      newcontig<-positive
    }
  }else{
    if(nrow(negative)>1){
      newcontig<-negative
    }
  }
    for(i in seq(1:nrow(newcontig))){
      tmp<-data[i,]
      data[i,]<-data[row<-which(data == as.character(newcontig["name"][i,])),]
      data[row,]<-tmp
    }
}
data["V5"]<-data["V4"]<-1
colnames(data)<- c("chr", "start", "end","V4","V5")
colnames(rnadata)<- c("Gene", "start", "end","Contig","color")
colnames(fdata)<-c("Gene","start","end","Contig","Cog")
colnames(rdata)<-c("Gene","start","end","Contig","Cog")
colnames(gcdata)<-c("Contig","start","end","gccontent","gcskew")
tocirmeasure<-segAnglePo(measure, seg=c(as.matrix(measure["seg.name"][,])))
tocirmeasure[2,2]<-300;tocirmeasure[2,3]<-620
tocirmeasure[3,2]<-0;tocirmeasure[3,3]<-620
tocirmeasure[4,2]<-420;tocirmeasure[4,3]<-620
tocirmeasure[5,2]<-90;tocirmeasure[5,3]<-620
tocirmeasure[6,2]<-120;tocirmeasure[6,3]<-620
tocirmeasure[7,2]<-180;tocirmeasure[7,3]<-620
tocirmeasure[8,2]<-240;tocirmeasure[8,3]<-620
tocirmeasure2<-segAnglePo(measure2, seg=c(as.matrix(measure2["seg.name"][,])))
tocirmeasure2[1,2]<-330
tocirmeasure4<-segAnglePo(measure4, seg=c(as.matrix(measure4["seg.name"][,])))
tocirmeasure4[1,2]<-395
tocirmeasure8<-segAnglePo(measure8, seg=c(as.matrix(measure8["seg.name"][,])))
tocirmeasure8[1,2]<-146
tocirmeasure10<-segAnglePo(measure10, seg=c(as.matrix(measure10["seg.name"][,])))
tocirmeasure10[1,2]<-210
tocir <- segAnglePo(data, seg=c(as.matrix(data["chr"][,])))
tocirf<-segAnglePo(fdata, seg=c(as.matrix(fdata["Gene"][,])))
tocirr<-segAnglePo(rdata, seg=c(as.matrix(rdata["Gene"][,])))
tocirrna<-segAnglePo(rnadata, seg=c(as.matrix(rnadata["Gene"][,])))
tocirgc<-segAnglePo(gcdata, seg=c(as.matrix(gcdata["Contig"][,])))
getl<-function (contig) {
  row<-which(data == as.character(contig))
  return(as.numeric(data["end"][row,]))
}
getal<-function (contig) {
  row<-which(tocir == as.character(contig))
  return(as.numeric(tocir[row,3]) -  as.numeric(tocir[row,2]))
}
getsa<-function (contig) {
  row<-which(tocir == as.character(contig))
  return(as.numeric(tocir[row,2]))
}
for(i in seq(1:nrow(tocirf))){
  tocirf[i,2]<-getsa(c(as.matrix(fdata["Contig"][i,])))+(getal(fdata["Contig"][i,])/getl(c(as.matrix(fdata["Contig"][i,]))))*as.numeric(fdata["start"][i,])
  tocirf[i,3]<-getsa(fdata["Contig"][i,])+(getal(fdata["Contig"][i,])/getl(c(as.matrix(fdata["Contig"][i,]))))*as.numeric(fdata["end"][i,])
  tocirf[i,1]<-paste(c("bar_"),i,sep="")
}
fdata["Gene"]<-as.data.frame(tocirf[,1])
for(i in seq(1:nrow(tocirr))){
  tocirr[i,2]<-getsa(rdata["Contig"][i,])+(getal(rdata["Contig"][i,])/getl(c(as.matrix(rdata["Contig"][i,]))))*as.numeric(rdata["start"][i,])
  tocirr[i,3]<-getsa(rdata["Contig"][i,])+(getal(rdata["Contig"][i,])/getl(c(as.matrix(rdata["Contig"][i,]))))*as.numeric(rdata["end"][i,])
  tocirr[i,1]<-paste(c("bar_"),i,sep="")
  }
rdata["Gene"]<-as.data.frame(tocirr[,1])
for(i in seq(1:nrow(tocirrna))){
  tocirrna[i,2]<-getsa(rnadata["Contig"][i,])+(getal(rnadata["Contig"][i,])/getl(c(as.matrix(rnadata["Contig"][i,]))))*as.numeric(rnadata["start"][i,])
  tocirrna[i,3]<-getsa(rnadata["Contig"][i,])+(getal(rnadata["Contig"][i,])/getl(c(as.matrix(rnadata["Contig"][i,]))))*as.numeric(rnadata["end"][i,])
  tocirrna[i,1]<-paste(c("bar_"),i,sep="")
}
rnadata["Gene"]<-as.data.frame(tocirrna[,1])
for(i in seq(1:nrow(tocirgc))){
  tocirgc[i,2]<-getsa(gcdata["Contig"][i,])+(getal(gcdata["Contig"][i,])/getl(c(as.matrix(gcdata["Contig"][i,]))))*as.numeric(gcdata["start"][i,])
  tocirgc[i,3]<-getsa(gcdata["Contig"][i,])+(getal(gcdata["Contig"][i,])/getl(c(as.matrix(gcdata["Contig"][i,]))))*as.numeric(gcdata["end"][i,])
  tocirgc[i,6]<-gcdata["start"][i,]
  tocirgc[i,7]<-gcdata["end"][i,]
  tocirgc[i,1]<-paste(c("bar_"),i,sep="")
  }
gcdata["Contig"]<-as.data.frame(tocirgc[,1])
if(makecog==c("Y")){
  colorlist<-list( "A" = "darkcyan", "B" = "sienna", "C" = "orange", "D" = "yellow", "E" = "green",
                   "F" = "pink", "G" = "steelblue", "H" = "purple", "I" = "brown", "J" = "violet", 
                   "K" = "maroon", "L" = "gold", "M" = "dark green", "N" = "dark red", "O" = "dark blue", 
                   "P" = "cyan", "Q" = "dark orange", "R" = "turquoise", "S" = "dark violet", "T" = "chocolate",
                   "U" = "beige", "V" = "yellowgreen", "W" = "dimgray", "X" = "gray" ,"Y" = "orangered", 
                   "Z" = "dark cyan")
  
  fdata["V6"]<-rdata["V6"]<-1
  fdata["V6"]<-as.character(colorlist[c(as.matrix(fdata["Cog"]))])
  rdata["V6"]<-as.character(colorlist[c(as.matrix(rdata["Cog"]))])
}
mymedian<-median(gcdata$gccontent)
pdf(file="Rplot.pdf", width = 10, height =10)
par(mar=c(2,2,2,2))
plot(c(0,1000), c(0,1000), type="n", axes=FALSE, xlab="", ylab="", main="")
if(nrow(data)>1){
  circos(R=400, cir=tocir, W=10,type="chr", print.chr.lab=F, scale=F)
}
if(makecog==c("Y")){
  circos(R=340, cir=tocirf, W=60, mapping=fdata, type="b3", col.v=7, col=c(as.matrix(fdata["V6"])), B=F,scale=F, lwd=abs(as.matrix((fdata["end"]-fdata["start"])/4000)))
  circos(R=290, cir=tocirr, W=60, mapping=rdata, type="b3", col.v=7, col=c(as.matrix(rdata["V6"])), B=F,scale=F, lwd=abs(as.matrix((rdata["end"]-rdata["start"])/4000)))
}else{
  circos(R=340, cir=tocirf, W=60, mapping=fdata, type="b3", col.v=7, col=c("dark green"), B=F,scale=F, lwd=abs(as.matrix((fdata["end"]-fdata["start"])/4000)))
  circos(R=290, cir=tocirr, W=60, mapping=rdata, type="b3", col.v=7, col=c("dark red"), B=F,scale=F, lwd=abs(as.matrix((rdata["end"]-rdata["start"])/4000)))
}
circos(R=250, cir=tocirrna, W=50, mapping=rnadata, type="b3", col.v=6, col=c(as.matrix(rnadata["color"])), B=F,scale=T, lwd=abs(as.matrix((rnadata["end"]-rnadata["start"])/4000)))
circos(R=200, cir=tocirgc, W=30, mapping=gcdata, type="b", col.v=4, col=ifelse(gcdata["gccontent"]>=mymedian,c("black"),c("white")), B=F,scale=F, lwd=0.4)
circos(R=200, cir=tocirgc, W=-30, mapping=gcdata, type="b", col.v=4, col=ifelse(gcdata["gccontent"]>=mymedian,c("white"),c("gray")), B=F,scale=F, lwd=0.4)
circos(R=130, cir=tocirgc, W=40, mapping=gcdata, type="b", col.v=5, col=ifelse(gcdata["gcskew"]>0,c("green"),c("white")), B=F,scale=F, lwd=0.4)
circos(R=130, cir=tocirgc, W=-40, mapping=gcdata, type="b", col.v=5, col=ifelse(gcdata["gcskew"]>0,c("white"),c("purple")), B=F,scale=F, lwd=0.4)
circos(R=62, cir=tocirmeasure, W=5, mapping=measure, type="b3", col.v=4, col=c("gray"), B=F, lwd=2)
circos(R=65, cir=tocirmeasure2, W=30, mapping=measure2, type="label2", col.v=4, col=c("gray"), B=F,scale=F, cex =0.5, side="out")
circos(R=65, cir=tocirmeasure4, W=30, mapping=measure4, type="label2", col.v=4, col=c("gray"), B=F,scale=F, cex =0.5, side="out")
circos(R=65, cir=tocirmeasure8, W=30, mapping=measure8, type="label2", col.v=4, col=c("gray"), B=F,scale=F, cex =0.5, side="out")
circos(R=65, cir=tocirmeasure10, W=30, mapping=measure10, type="label2", col.v=4, col=c("gray"), B=F,scale=F, cex =0.5, side="out")
legend(350,420, text.col=c("gray"), cex=0.6,legend="x Kb", box.col="white") ;
if(makecog== c("Y")){
  code<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
  cogletters<-rbind(fdata["Cog"],rdata["Cog"])
  letters<-ifelse(code %in% c(as.matrix(cogletters)),code,FALSE)
  letters<-letters[which(letters != FALSE)]
  lcolors<-as.character(colorlist[c(as.matrix(letters))])
  
  legend("right", legend = c(letters, "tRNA","rRNA", "GC Content","GC Skew +","GC Skew -"), 
         ncol = 1,
         xpd = NA, cex = 0.8,  bty="n",
         fill=c(lcolors,"red","blue","black","purple","green"),
         border = c("white"),
         title = "COG Categories") 
}else{
  legend("right", legend = c("Forward genes","Reverse Genes","tRNA","rRNA", "GC Content","GC Skew +","GC Skew -"), 
         ncol = 1,
         xpd = NA, cex = 0.8,  bty="n",
         fill=c("dark green","dark red","red","blue","black","purple","green"),
         border = c("white"))  
}
dev.off()
""")

	plotstep.close


def MakeCog(gbk):
	wd=os.getcwd()
	gbkname=gbk.replace("/"," ").split()[len(gbk.replace("/"," ").split())-1]
	protdict={}
	contiggene={}
	recs = [rec for rec in SeqIO.parse(gbk, "genbank")]
	for rec in recs:
		contigname=rec[0:].id
		feats = [feat for feat in rec.features if feat.type == "CDS"]

		for feat in feats:
			cdsname=str(feat.qualifiers["locus_tag"]).replace("'","").replace("[","").replace("]","")
			if "translation" in feat.qualifiers:	
				translation=str(feat.qualifiers["translation"]).replace("'","").replace("[","").replace("]","")
				band=0
				for seq in protdict.items():
					if seq[1] == translation:
						band=1

				if band==0:
					protdict[cdsname]=translation
					contiggene[cdsname]=contigname

	faa= open(str(wd+"/"+gbkname+ ".faa"), 'w')
	for cds in protdict.items():
		faa.write(">%s | %s\n" % (cds[0], contiggene[cds[0]]))
		faa.write("%s\n" % (cds[1]))
	
	faa.close
	return str(gbkname + ".faa")

def callRPSBlast(faafile,threads):
	faaname=faafile.replace("/"," ").split()[len(faafile.replace("/"," ").split())-1]
	if os.path.isfile(str("rpsblast."+faaname+".out")) == False:
		subprocess.call(["data/rpsblast", "-query", faafile, "-db", "data/Cog_LE/Cog", "-out", str("rpsblast."+faaname+".out"), "-evalue", "1e-2", "-outfmt", "6", "-num_threads",str(threads)])
	
	subprocess.call(["rm","-rf",str("results_"+faaname)])
	subprocess.call(["perl", "data/cdd2cog.pl", "-r", str("rpsblast."+faaname+".out"), "-c","data/cddid.tbl","-f", "data/fun.txt" ,"-w", "data/whog", "-a"])
	subprocess.call(["mv","results",str("results_"+faaname)])
	subprocess.call(["mv",str(faaname),str("results_"+faaname)])

	os.chdir(str("results_"+faaname))
	parsedcog= open("parsedcog.dat", 'w')
	with open("rps-blast_cog.txt") as tsv:
		rpsblastout=csv.reader(tsv, delimiter="\t")
		first_row = next(rpsblastout)
		first_row = next(rpsblastout)
		prevname=first_row[0]
		preveval=first_row[10]
		prevcog=first_row[13]
		for line in rpsblastout:

			if prevname!=line[0]:
				parsedcog.write(str(prevname+" "+prevcog+"\n"))
				prevname=line[0]
				preveval=line[10]
			else:
				if preveval>line[10]:
					prevname=line[0]
					preveval=line[10]
					prevcog=line[13]


	parsedcog.write(str(prevname+" "+prevcog+"\n"))
	parsedcog.close
	#return str("results_"+faaname)

def getCog(genename):
		
		cog=""
		toplot=open("parsedcog.dat",'r')
		for line in toplot:
			if re.search(genename, line):
				cog=line.split()[1]  #0 gene name, 1 cog letter
		
		if cog == "":
			cog="S"

		toplot.close
		return cog

def GBKParser(genbank_file, makecog, filterc):
	gbkname=genbank_file.replace("/"," ").split()[len(genbank_file.replace("/"," ").split())-1]
	fna=open(str(gbkname+".fna"), 'w')
	contig=open("contigplot.dat", 'w')#only contgis that contains genes
	forward = open("forwardplot.dat", 'w')
	reverse = open("reverseplot.dat", 'w')
	rna= open("rna.dat", 'w')

	genomelength=0
	recs = [rec for rec in SeqIO.parse(str("../"+genbank_file), "genbank")]  
	for rec in recs:
		contigname=rec[0:].id
		contiglen=len(rec[0:].seq)
		contigseq=rec[0:].seq
		genomelength+=contiglen
		contigband=0
		feats = [feat for feat in rec.features if feat.type == "CDS"]
		for feat in feats:
			feat.location=str(feat.location).replace("("," ").replace(")","").replace("]","").replace("[","").replace(">","").replace("<","").replace(":"," ")
			locations=feat.location.split()
			gene=str(feat.qualifiers["locus_tag"]).replace("'","").replace("[","").replace("]","")

			if makecog=="Y":
				cog=getCog(gene)
				contigband=1
				if locations[2]=="+":
					forward.write("%s %s %s %s %s" % (gene, int(locations[0])+1, int(locations[1])+1, contigname, cog))  #gene, start, end, contig
					forward.write("\n")
				else:
					reverse.write("%s %s %s %s %s" % (gene, int(locations[0])+1, int(locations[1])+1, contigname, cog))  #gene, start, end, contig
					reverse.write("\n")
			else:
				contigband=1
				if locations[2]=="+":
					forward.write("%s %s %s %s %s" % (gene, int(locations[0])+1, int(locations[1])+1, contigname, "E"))  #just for color
					forward.write("\n")
				else:
					reverse.write("%s %s %s %s %s" % (gene, int(locations[0])+1, int(locations[1])+1, contigname, "M"))  #just for color
					reverse.write("\n")

		rrnafeats = [feat for feat in rec.features if feat.type == "rRNA"]
		for feat in rrnafeats:
			feat.location=str(feat.location).replace("("," ").replace(")","").replace("]","").replace("[","").replace(">","").replace("<","").replace(":"," ")
			locations=feat.location.split()
			gene=str(feat.qualifiers["locus_tag"]).replace("'","").replace("[","").replace("]","")
			

			contigband=1
			rna.write("%s %s %s %s %s" % (gene, int(locations[0])+1, int(locations[1])+1, contigname, "blue"))  #gene, start, end, contig
			rna.write("\n")

		trnafeats = [feat for feat in rec.features if feat.type == "tRNA"]
		for feat in trnafeats:
			feat.location=str(feat.location).replace("("," ").replace(")","").replace("]","").replace("[","").replace(">","").replace("<","").replace(":"," ")
			locations=feat.location.split()
			gene=str(feat.qualifiers["locus_tag"]).replace("'","").replace("[","").replace("]","")
			
			contigband=1
			rna.write("%s %s %s %s %s" % (gene, int(locations[0])+1, int(locations[1])+1, contigname, "red"))  #gene, start, end, contig
			rna.write("\n")

		if filterc == "Y":
			if contigband==1:
				contig.write(str(contigname+" "+"1 "+str(contiglen)))
				contig.write("\n")
		else:
			contig.write(str(contigname+" "+"1 "+str(contiglen)))
			contig.write("\n")

		fna.write(">%s \n" % (contigname))
		fna.write("%s\n" % (contigseq))




	fna.close
	contig.close
	forward.close
	reverse.close
	rna.close
	return genomelength

def GC_content_window(s):
	gc = sum(s.count(x) for x in ['G','C','g','c','S','s'])
	gc_content = gc/float(len(s))
	return round(gc_content,4) 


def GC_skew_window(s):
	g = s.count('G')+s.count('g')
	c = s.count('C')+s.count('c')

	try:
		skew = (g-c)/float(g+c)
	except ZeroDivisionError:
		skew = 0
	return round(skew,4)


def GCcalc(filename,window,step,filteredcontigs):

	filterc=open(filteredcontigs,'r')
	filteredc=filterc.read()
	gcplot=open("gcskewplot.dat", 'w')#only contgis that contains genes

	seqobj = SeqIO.parse(filename,'fasta')

	for record in seqobj: 
		pos_array = []
		gc_content_value_array = []
		gc_skew_value_array = []
		name = record.id
		seq = record.seq
		start = 0
		end = 0
		gc = 0
		gc_skew = 0
		if name in filteredc:
			for i in range(0,len(seq),step):
				subseq = seq[i:i+window]
				gc_content = (GC_content_window(subseq))
				gc_skew = (GC_skew_window(subseq))
				start = (i + 1 if (i+1<=len(seq)) else i)
				end = ( i + step if (i+ step<=len(seq)) else len(seq))
				gcplot.write("%s %s %s %s %s\n" % (name,start,end,gc_content,gc_skew))

		gcplot.close
		filterc.close

def main():
	parser = OptionParser(usage = "Usage: python wrapper.py -f genbankfile.gbk")
	parser.add_option("-C","--CogAssign",dest="cogoption",help="default:N, Y or N assign cog function by rpsblast", default="N")
	parser.add_option("-F","--FilterContigs",dest="filterc",help="default:Y, show only contigs that contains genes",default="Y")
	parser.add_option("-f","--file",dest="filename",help="Input Fasta format file",metavar="GENBANK FILE")
	parser.add_option("-w","--window",dest="window",help="default:3000, window to take for gccontent and gc skew",default=3000)
	parser.add_option("-s","--step",dest="step",help="default:1500 step to move your window",default=1500)
	parser.add_option("-t","--threads",dest="threads",help="default:1 threads for rpsblast",default=1)


	(options,args) = parser.parse_args()

	makecog = options.cogoption
	filterc= options.filterc
	genbank_file = options.filename
	window = int(options.window)
	step = int(options.step)
	threads = int(options.threads)
	gbkname=genbank_file.replace("/"," ").split()[len(genbank_file.replace("/"," ").split())-1]


	if makecog == "Y":
		print "Making cog assign"
		faafile=MakeCog(genbank_file)
		callRPSBlast(faafile,threads)
	else:
		resultsfolder=str("results_"+gbkname+".faa")
		subprocess.call(["rm","-rf",str("results_"+gbkname+".faa")])
		subprocess.call(["mkdir",str("results_"+gbkname+".faa")])
		os.chdir(resultsfolder)

	print "Parsing GBK"
	gl=GBKParser(genbank_file, makecog, filterc)
	print "Computing GC content and GC Skew +/-"
	GCcalc(str(gbkname+".fna"),window,step,"contigplot.dat")
	print "Plotting"
	printPlotStep()
	subprocess.call(["Rscript", "plotstep.R", str(gl), str(makecog)])
	subprocess.call(["mv","Rplot.pdf",str(str(gbkname).replace(".gbk","")+".pdf")])
	subprocess.call(["rm","-rf","plotstep.R"])

	print "Done :D"



if __name__ == '__main__':
	main()
	sys.exit()