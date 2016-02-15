rm(list=ls());
library(Rsamtools)
library(OmicCircos)
library(data.table)
args<-commandArgs()
gl<-as.numeric(args[6])
makecog<-c(args[7])
#gl<-5827046
#setwd("results_putida.gbk.faa/")
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

#ORDER GCSKEW DATA
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


#head(data)
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
#tocirgc<-segAnglePo(gcdata, seg=c(as.matrix(gcdata["Gene"][,])))
tocirgc<-segAnglePo(gcdata, seg=c(as.matrix(gcdata["Contig"][,])))
#head(tocirgc)

getl<-function (contig) {
  row<-which(data == as.character(contig))
  return(as.numeric(data["end"][row,]))
}

#get angle length of contig
getal<-function (contig) {
  row<-which(tocir == as.character(contig))
  return(as.numeric(tocir[row,3]) -  as.numeric(tocir[row,2]))
}

getsa<-function (contig) {
  row<-which(tocir == as.character(contig))
  return(as.numeric(tocir[row,2]))
}

#adjust angle genes for each contig, if you get an error, maybe some contig not exist in data
#forward genes
for(i in seq(1:nrow(tocirf))){
  tocirf[i,2]<-getsa(c(as.matrix(fdata["Contig"][i,])))+(getal(fdata["Contig"][i,])/getl(c(as.matrix(fdata["Contig"][i,]))))*as.numeric(fdata["start"][i,])
  tocirf[i,3]<-getsa(fdata["Contig"][i,])+(getal(fdata["Contig"][i,])/getl(c(as.matrix(fdata["Contig"][i,]))))*as.numeric(fdata["end"][i,])
  tocirf[i,1]<-paste(c("bar_"),i,sep="")
  
}
fdata["Gene"]<-as.data.frame(tocirf[,1])

#reverse genes
for(i in seq(1:nrow(tocirr))){
  tocirr[i,2]<-getsa(rdata["Contig"][i,])+(getal(rdata["Contig"][i,])/getl(c(as.matrix(rdata["Contig"][i,]))))*as.numeric(rdata["start"][i,])
  tocirr[i,3]<-getsa(rdata["Contig"][i,])+(getal(rdata["Contig"][i,])/getl(c(as.matrix(rdata["Contig"][i,]))))*as.numeric(rdata["end"][i,])
  tocirr[i,1]<-paste(c("bar_"),i,sep="")
  
  }
rdata["Gene"]<-as.data.frame(tocirr[,1])

#rna genes (trna rrna)
for(i in seq(1:nrow(tocirrna))){
  tocirrna[i,2]<-getsa(rnadata["Contig"][i,])+(getal(rnadata["Contig"][i,])/getl(c(as.matrix(rnadata["Contig"][i,]))))*as.numeric(rnadata["start"][i,])
  tocirrna[i,3]<-getsa(rnadata["Contig"][i,])+(getal(rnadata["Contig"][i,])/getl(c(as.matrix(rnadata["Contig"][i,]))))*as.numeric(rnadata["end"][i,])
  tocirrna[i,1]<-paste(c("bar_"),i,sep="")
  
}
rnadata["Gene"]<-as.data.frame(tocirrna[,1])

#gc
for(i in seq(1:nrow(tocirgc))){
  tocirgc[i,2]<-getsa(gcdata["Contig"][i,])+(getal(gcdata["Contig"][i,])/getl(c(as.matrix(gcdata["Contig"][i,]))))*as.numeric(gcdata["start"][i,])
  tocirgc[i,3]<-getsa(gcdata["Contig"][i,])+(getal(gcdata["Contig"][i,])/getl(c(as.matrix(gcdata["Contig"][i,]))))*as.numeric(gcdata["end"][i,])
  tocirgc[i,6]<-gcdata["start"][i,]
  tocirgc[i,7]<-gcdata["end"][i,]
  tocirgc[i,1]<-paste(c("bar_"),i,sep="")
  }
gcdata["Contig"]<-as.data.frame(tocirgc[,1])

if(makecog==c("Y")){
  getcolor<-function(type){
    switch(as.character(type),
           "A" = "darkcyan", "B" = "sienna", "C" = "orange", "D" = "yellow", "E" = "green", 
           "F" = "pink", "G" = "steelblue", "H" = "purple", "I" = "brown", "J" = "violet", 
           "K" = "maroon", "L" = "gold", "M" = "dark green", "N" = "dark red", "O" = "dark blue", 
           "P" = "cyan", "Q" = "dark orange", "R" = "turquoise", "S" = "gray", "T" = "chocolate",
           "U" = "beige", "V" = "yellowgreen", "W" = "dimgray", "Y" = "orangered", 
           "Z" = "dark cyan")
  }
  fdata["V6"]<-rdata["V6"]<-1
  for(i in seq(1:nrow(as.data.frame(fdata)))){
    fdata["V6"][i,]<-getcolor(fdata["Cog"][i,])
  }
  for(i in seq(1:nrow(as.data.frame(rdata)))){
    rdata["V6"][i,]<-getcolor(rdata["Cog"][i,])
  }
}


#to gcskew
#idx = pmatch(c(as.matrix(gcdata["Gene"])), names(gr))  ## NAs if duplicates or not unique
#seq = getSeq(fa, gr[as.vector(as.matrix(idx))])
#freq<-alphabetFrequency(seq, baseOnly=TRUE)
#avggc<-sum((freq[,2]+freq[,3])/(freq[,1]+freq[,2]+freq[,3]+freq[,4]))/nrow(freq)
#gccol<-data.frame(gccontent=seq(1:nrow(as.data.frame(seq))))

#gccol["gccontent"]<-(freq[,2]+freq[,3])/(freq[,1]+freq[,2]+freq[,3]+freq[,4])


#getting fluctuation
#for(i in seq(1:nrow(as.data.frame(seq)))){
#  gccol["gccontent"][i,]<-c(avggc-((freq[i,2]+freq[i,3])/((freq[i,1]+freq[i,2]+freq[i,3]+freq[i,4]))))
#}
#gcdata["gccontent"]<-gccol
######gc skew forward+reverse ###########################
#gcdata["V6"]<-1
#freq<-alphabetFrequency(seq, baseOnly=TRUE)
#gcdata["V6"]<-(freq[,3]-freq[,2])/(freq[,3]+freq[,2])

#########################################################
#####intento de lineas

#head(gcdata);head(gccontent)
#gccontent<-data.frame(seg.name=seq(1:(nrow(gcdata)*2)), seg.po=1, gccontent=1, gcskew=1)
#j=1
#contig<-c("foo")
#for(i in seq(1:nrow(gcdata))){
#  gccontent["seg.name"][j,]<-c(as.matrix(gcdata["Gene"][i,]))
#  gccontent["seg.name"][j+1,]<-c(as.matrix(gcdata["Gene"][i,]))
#  gccontent["seg.po"][j,]<-c(as.matrix(gcdata["start"][i,]))
#  gccontent["seg.po"][j+1,]<-c(as.matrix(gcdata["end"][i,]))
  
#  if(contig!=c(as.matrix(gcdata["Contig"][i,]))){
#    contig<-c(as.matrix(gcdata["Contig"][i,]))
#    gccontent["gccontent"][j,]<-0
#    gccontent["gccontent"][j+1,]<-gcdata["gccontent"][i,]
    #gccontent["gcskew"][j,]<-0
    #gccontent["gcskew"][j+1,]<-c(as.matrix(gcdata["V6"][i,]))

#  }else{
#    gccontent["gccontent"][j,]<-gcdata["gccontent"][i-1,]
#    gccontent["gccontent"][j+1,]<-gcdata["gccontent"][i,]
    #gccontent["gcskew"][j,]<-c(as.matrix(gcdata["V6"][i-1,]))
    #gccontent["gcskew"][j+1,]<-c(as.matrix(gcdata["V6"][i,]))
#  }
#  j<-j+2
#}

##########
#range(gcdata["gccontent"])
#sum(gcdata["gccontent"])/nrow(gcdata)
#plot(density(gcdata$gccontent))
mymedian<-median(gcdata$gccontent)
#mean(gcdata$gccontent)
#class(gcdata["gccontent"])
#head(gcdata)
#range(gcdata["gccontent"])
#range(gcdata["gcskew"])

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
legend(370,420, text.col=c("gray"), cex=0.6,legend="x Kb", box.col="white") ;


if(makecog== c("Y")){
  legend("right", legend = c("A ","B ", "C ", "D ", "E ", "F ","G ","H ", "I ", "J ", "K ","L ","M ", "N ","O ", "P ","Q ",
                             "R ","T ", "U ","V ","W ", "Y ","Z ","Unknown","tRNA","rRNA", "GC Content","GC Skew +","GC Skew -"), 
         ncol = 1,
         xpd = NA, cex = 0.75,  bty="n",
         fill=c("darkcyan","sienna", "orange","yellow","green","pink","steelblue","purple","brown","violet",
                "maroon","gold","dark green","dark red","dark blue","cyan","dark orange","turquoise",
                "chocolate","beige","yellowgreen","dimgray","orangered","dark cyan","gray","red","blue","black","purple","green"),
         border = c("white"),
         title = "COG Categories") 
}else{
  legend("right", legend = c("Forward genes","Reverse Genes","tRNA","rRNA", "GC Content","GC Skew +","GC Skew -"), 
         ncol = 1,
         xpd = NA, cex = 0.75,  bty="n",
         fill=c("dark green","dark red","red","blue","black","purple","green"),
         border = c("white"))  
}

dev.off()

