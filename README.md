![banner](https://raw.githubusercontent.com/microgenomics/tutorials/master/img/microgenomics.png)
# plotMyGBK
-----------
plotMyGBK is a very simple to use pipe that use a genbank file and give you a nice circular graph (using omiccircos to make the graph)

## Requisites

* Python >= 2.7 with the module [Biopython](http://biopython.org/wiki/Download)
* R (tested in 3.2.3), and the following libraries: 
	* [Rsamtool](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html) 
	* [OmicCircos](http://bioconductor.org/packages/release/bioc/html/OmicCircos.html)
	* [data.table](https://cran.r-project.org/web/packages/data.table/index.html)

##Installation

* If you want to classify your genes by functions, download the zip from [here](https://www.dropbox.com/s/n6ycyhrtlz3ixeb/data.zip?dl=0), and unzip. This package have all neccesary files to classify your genes using the Cog Assignments.
* Download the script from [this repository](=) and put them together with data folder that you unzipped previously.

## Usage

plotMyGBK have minimal use:

	python plotMyGBK.py -f mygenbank.gbk

and a complete use:
	
	python plotMyGBK.py -f mygenbank.gbk -C Y -F Y -w 3000 -s 1500

* -f is the genbank file
* -C Cog Assign, this will call to rpsblast (provided in data folder, v=2.2.31+), and classify the genes by Cog function. default=N
* -F Filter Contigs, show only contigs that have genes. default=Y
* -w Windows size, to computing GC content and GC skew set the windows size, default:3000
* -s Step size, size in bp to move windows (-w), default:1500

##Trick

plotMyGBK first search for rpsblast out (format 6), if this file exists, you can save a lot of time!, put this file together with the scripts and change the name to: rpsblast.GBKNAME.faa.out where GBKNAME is the name of your gbk file included the extension .gbk (e.g. mygbk.gbk ==> rpsblastout.mygbk.gbk.faa.out)

##Output

Since your gbk, plotMyGBK will generate a folder called results\_GBKNAME.faa where GBKNAME is the genbank file, this folder contain several files like faa, fna, parsed rpsblast results, other files necessaries to plot and last the plot in pdf format.

##Notes
* When you execute plotMyGBK.py make sure data folder and plotstep.R are in that directory.
* The rpsblast in the data folder was compiled in Mac OSX, if you use other OS, you have to put the corresponding binary of rpsblast in this data folder