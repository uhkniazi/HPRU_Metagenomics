HPRU_Metagenomics
=================
  
Pipeline for metagenomics data analysis.

# Requirements
## 1- Diamond Blast
[From the diamond blast webpage](http://ab.inf.uni-tuebingen.de/software/diamond/) DIAMOND is a new alignment tool for aligning short DNA sequencing reads to a protein reference database such as NCBI-NR. On Illumina reads of length 100-150bp, in fast mode, DIAMOND is about 20,000 times faster than BLASTX, while reporting about 80-90% of all matches that ...
  
[Tutorial](http://ab.inf.uni-tuebingen.de/data/software/diamond/download/public/manual.pdf) Opens a PDF document.
  
#### Download and Installation
This assumes you have a linux bash shell open.

```
wget http://github.com/bbuchfink/diamond/releases/download/v0.8.8/diamond-linux64.tar.gz
tar -xvf diamond-linux64.tar.gz
```
#### Preparing the blast nr database
Around 60 GB of more of diskspace may be required.
```
nohup wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz &
nohup ./diamond makedb --in /Capmeta/home/Downloads/blastdb/nr.fa -d nr &
```
  
## 2- MEGAN
[From the MEGAN webpage](http://ab.inf.uni-tuebingen.de/software/megan6/)MEGAN6 is a comprehensive toolbox for interactively analyzing microbiome data. All the interactive tools you need in one ...

**Please Note** this pipeline was tested with MEGAN5, so you can apply the same steps to MEGAN6, assuming it all works.
#### Download and Installation
```
wget http://ab.inf.uni-tuebingen.de/data/software/megan5/download/MEGAN_unix_5_11_3.sh
chmod +x MEGAN_unix_5_11_3.sh
./MEGAN_unix_5_11_3.sh
```
**Please Note** [From MEGAN5 Webpage](http://ab.inf.uni-tuebingen.de/data/software/megan5/download/welcome.html) Auxiliary mapping files for taxon analysis:

If the BLAST or other comparison files do not contain taxon names, but do contain GI numbers, then load the following mapping file to allow MEGAN to determine taxon names. The file has to be unzipped before it can be used:

Download the database using
```
wget http://www-ab.informatik.uni-tuebingen.de/data/software/megan5/download/gi_taxid-March2015X.zip
unzip gi_taxid-March2015X.zip
```
