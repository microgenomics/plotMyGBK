import sys, os, re, multiprocessing
import subprocess, csv
from optparse import OptionParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO


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

def callRPSBlast(faafile):
	faaname=faafile.replace("/"," ").split()[len(faafile.replace("/"," ").split())-1]
	if os.path.isfile(str("rpsblast."+faaname+".out")) == False:
		subprocess.call(["data/rpsblast", "-query", faafile, "-db", "data/Cog_LE/Cog", "-out", str("rpsblast."+faaname+".out"), "-evalue", "1e-2", "-outfmt", "6", "-num_threads",str(multiprocessing.cpu_count())])
	
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
				if locations[2]=="+":
					forward.write("%s %s %s %s" % (gene, int(locations[0])+1, int(locations[1])+1, contigname))  #gene, start, end, contig
					forward.write("\n")
				else:
					reverse.write("%s %s %s %s" % (gene, int(locations[0])+1, int(locations[1])+1, contigname))  #gene, start, end, contig
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

	(options,args) = parser.parse_args()

	makecog = options.cogoption
	filterc= options.filterc
	genbank_file = options.filename
	window = options.window
	step = options.step
	gbkname=genbank_file.replace("/"," ").split()[len(genbank_file.replace("/"," ").split())-1]


	if makecog == "Y":
		print "Making cog assign"
		faafile=MakeCog(genbank_file)
		callRPSBlast(faafile)
	else:
		subprocess.call(["mkdir",str("results_"+genbank_file)])
		resultsfolder=str("results_"+gbkname)
		os.chdir(resultsfolder)

	print "Parsing GBK"
	gl=GBKParser(genbank_file, makecog, filterc)
	print "Computing GC content and GC Skew +/-"
	GCcalc(str(gbkname+".fna"),window,step,"contigplot.dat")
	print "Plotting"
	subprocess.call(["Rscript", "../plotstep.R", str(gl), str(makecog)])
	subprocess.call(["mv","Rplot.pdf",str(str(gbkname).replace(".gbk","")+".pdf")])

	print "Done :D"



if __name__ == '__main__':
	main()
	sys.exit()