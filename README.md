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
[From the MEGAN webpage](http://ab.inf.uni-tuebingen.de/software/megan6/) MEGAN6 is a comprehensive toolbox for interactively analyzing microbiome data. All the interactive tools you need in one ...

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
# Running the pipeline
The starting point is a FASTQ file from the sequencing run while the output
file at the end will be a count table produced by MEGAN.

*FASTQ File* --> **Diamond Blast** --> *.daa file* --> **Diamond Blast**
--> *.m8 blast format file* --> **MEGAN** --> *.rma file* --> **MEGAN**
--> *.txt count table* --> **Statistical Analysis**

## Running Diamond Blast
This is the most time consuming task and resource hungry task in the pipeline.
The script below is an example of how to run two samples through Diamond Blast.
```
# script to run diamond blast on several files
# these are the names of the .fastq files in the input directory
BM_samples="003_S1_L001_R1_001 003_S1_L001_R2_001"

for s in ${BM_samples}
do
    ./diamond blastx -d nr -q input/${s}.fastq -a BM_${s} -t temp/
done
```
Once the jobs are finished, respective *.daa* files will be produced in the current directory. Run the following script to produce blast tab separated .m8 format files

```
# script to run diamond blast on several files produced from previous script
BM_samples="BM_003_S1_L001_R1_001 BM_003_S1_L001_R2_001"

for s in ${BM_samples}
do
    ./diamond view -a ${s}.daa -o ${s}.m8
done
```
The various columns of the .m8 file stand for:

```
'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
qseqid means Query Seq-id
sseqid means Subject Seq-id
pident means Percentage of identical matches
length means Alignment length
mismatch means Number of mismatches
gapopen means Number of gap openings
qstart means Start of alignment in query
qend means End of alignment in query
sstart means Start of alignment in subject
send means End of alignment in subject
evalue means Expect value
bitscore means Bit score
```
## Running MEGAN5
MEGAN5 has to be run from the commandline, and a license also needs to be available in the same directory. Once a license has been obtained, copy the text
in the a text file e.g. License.txt.

**NOTE** MEGAN5 will run on one input file at a time, if two input files are provided, then it will merge the counts together from the two files
and produces one output file.

Create a file *commands.txt* in the same directory before running MEGAN5, a sample is shown below:
```
load taxGIFile='/Capmeta/home/Programs/megan/gi_taxid-March2015X.bin';
import blastFile='diamond/BM_003_S1_L001_R1_001.m8' meganFile='output/BM_003_S1_L001_R1_001.rma' maxMatches=100 minScore=60.0 maxExpected=1.0E-6 topPercent=10.0 minSupportPercent=0.01 minSupport=1 minComplexity=0.0 useMinimalCoverageHeuristic=false useSeed=false useCOG=false useKegg=false paired=false useIdentityFilter=false textStoragePolicy=Embed blastFormat=BlastTAB mapping='Taxonomy:GI_MAP=true';
quit;

```
Once the output file is produced, the count table can be produced using a second command file. **However for a cleaner picture of organism abundance it is recommended to cut the tree at *species level* and export the *summarized counts*.**  e.g. *export_commands_species.txt*:  
  
```
open file='output/S4.rma';
collapse rank='Species';
select rank='Species';
export what=CSV format=taxonname_count separator=comma counts=summarized file='output/S4.rma.txt'
;
quit;
```

If all the counts from the **root node** are required then use the following script e.g. *export_commands.txt*:
```
open file='output/BM_003_S1_L001_R1_001.rma';
unCollapse nodes=all;
select nodes=all;
export what=CSV format=taxonname_count separator=comma counts=assigned file='output/BM_003_S1_L001_R1_001.rma.txt';
quit;

```
MEGAN5 can be run from the command line using the following syntax, once the command files have been created:
```
xvfb-run --auto-servernum --server-num=1 ./MEGAN -g -E -L License.txt -c commands.txt
xvfb-run --auto-servernum --server-num=1 ./MEGAN -g -E -L License.txt -c export_commands.txt
```
In case we want to run two samples at a time, read the note above, the command file example is given below:
```
load taxGIFile='/Capmeta/home/Programs/megan/gi_taxid-March2015X.bin';
import blastFile='diamond/BM_004_S2_L001_R1_001.m8', 'diamond/BM_004_S2_L001_R2_001.m8' meganFile
='output/Sample_004_S2-nofa.rma' maxMatches=100 minScore=50.0 maxExpected=0.01 topPercent=10.0 mi
nSupport=1 minComplexity=0.0 useMinimalCoverageHeuristic=false useSeed=false useCOG=false useKegg
=false paired=false useIdentityFilter=false textStoragePolicy=Embed blastFormat=BlastTAB mapping=
'Taxonomy:GI_MAP=true';
quit;
```
### Changing MEGAN parameters
The input parameters for the *commands.txt, export_commands.txt* can be changed. The simplest way to know which command to use, is to run MEGAN in GUI mode and note down the commands that are displayed
in the terminal - then copy and paste these commands in the respective command file.

# Results Analysis
Once the csv count file is produced by MEGAN, this file can be analyzed in **R**:  

1- *megan_multi_plots.R*  
2- *megan_csv_import.R*  

**Note** both these script require the R package LearnBayes, to simulate from the dirichlet distribution. It should be installed within R using  
```R
install.packages("LearnBayes")
```

## megan_multi_plots.R
This script will run from the commandline and take 2 column csv file, produced by MEGAN using the script *export_commands_species.txt* or *export_commands.txt*.
It can take multiple input files, and will produce 5 output png plots in the same folder as the input file/s.  

```
Rscript megan_multi_plots.R Data_external2/S3.rma.txt Data_external2/S4.rma.txt
```
  
The above command will produce 2 sets of 5 plots, one for S3.rma.txt and second for S4.rma.txt. 


![alt text](https://github.com/uhkniazi/HPRU_Metagenomics/blob/master/Images/S3.rma.txt1.png "Top 10% of Abundant Taxa")  
**Figure 1** shows the proportion of reads/counts assigned by MEGAN to the top 10% taxa.  

![alt text](https://github.com/uhkniazi/HPRU_Metagenomics/blob/master/Images/S3.rma.txt2.png "Bar plot of top 10% Abundant Taxa")  
**Figure 2** shows the proportion of reads/counts assigned by MEGAN to the top 10% taxa in a Bar plot.  

![alt text](https://github.com/uhkniazi/HPRU_Metagenomics/blob/master/Images/S3.rma.txt3.png "Pie chart of top 10% Abundant Taxa")  
**Figure 3** shows the proportion of reads/counts assigned by MEGAN to the top 10% taxa in a Pie chart.   

![alt text](https://github.com/uhkniazi/HPRU_Metagenomics/blob/master/Images/S3.rma.txt4.png "Bar plot of top 30% Abundant Taxa")  
**Figure 4** shows the proportion of reads/counts assigned by MEGAN to the top 30% taxa in a Bar plot.  
  

**Figure 5** (not shown) is a Pie chart representation of the same data from Figure 4.


## megan_csv_import.R
This script will is run interactively, and is very similar to the previous script. The details of individual functions used in both scripts is described below:


### Functions used by script.
#### getAlpha
**ARGS**

**df** data.frame with 2 columns where column 2 (V2) contains the count data, and V1 the name of the species/taxa


**prior=1/2** Jeffery's non-informative prior for dirichlet distribution.

The function will remove the root node and return the alpha parameter for the dirichlet distribution.

**RETS** named vector

#### getPosterior
**ARGS**

**alpha** value of alpha for dirichlet distribution.

**n=1000** number of random samples to simulate from dirichlet distribution.

The function will simulate **n** samples from dirichlet posterior, using the value of **alpha** generated from **getAlpha** function.

**RETS** matrix

#### plot.bar

Simple function to plot a bar plot with 95% confidence interval error bars.

### Script Logic
1- Load the 2 column csv file produced by MEGAN.

2- Get the alpha for the posterior dirichlet, by using a non-informative dirichlet prior and a multinomial sampling model.  
3- Removed the group 'Not assgined' *feel free to comment out this line**.


```R
## data import
dfData = read.csv(file.choose(), stringsAsFactors = F, header=F)

ivAlpha = getAlpha(dfData)
# removed unassigned
# i = which(names(ivAlpha) == 'Not assigned')
# ivAlpha = ivAlpha[-i]
```
4- Simulate sample from dirichlet posterior to calculate confidence intervals.  
5- Get the mean (average) proportion for each taxa.  
6- Group this vector into 10 groups based on quantiles. Group 10 is the top 10% quantile, with perhaps the most interesting taxa.


```R
iAve = colMeans(mDir.post)
head(iAve)

# break into groups
groups = cut(iAve, breaks = quantile(iAve, 0:10/10), include.lowest = T, labels = 1:10)

plot(sort(iAve[groups == '10']), type='l', main='Top 10% Abundant Taxa',
     xlab='Taxa', ylab='Proportion')
```
7- Make a bar plot of the top 10% most abundant taxa.  
```R
cvTop = names(sort(iAve[groups == '10'], decreasing = T))

## plot the top 10% samples
mPlot = mDir.post[,cvTop]
plot.bar(mPlot, title = 'Sample 004_S1')
```

------------------------
# Metaphlan2

MetaPhlAn2 is a tool for profiling the composition of microbial communities from metagenomic shotgun sequencing data. In our tests, this tool is very fast for a quick analysis for the sample, however a more detailed approach using the pipeline mentioned earlier is recommended.


[Metaphlan2 tutorial](https://bitbucket.org/biobakery/metaphlan2)

## Installation
[Install](https://bitbucket.org/biobakery/metaphlan2#markdown-header-installation)

## Running
[Running](https://bitbucket.org/biobakery/metaphlan2#markdown-header-basic-usage)
```
# script to run metaphlan2 pipeline
# set environment variables
export PATH=`pwd`:$PATH
export mpa_dir=`pwd`

# two fastq files from samples
BM_samples="1467126010_clean 1467126010_cont"

for s in ${BM_samples}
do
    ./metaphlan2.py input/${s}.fq --bt2_ps very-sensitive-local --input_type multifastq --bowtie2out BM_${s}.bt2out > profiled_samples/BM_${s}.txt
done

# merge results tables the tables
utils/merge_metaphlan_tables.py profiled_samples/*.txt > output/merged_abundance_table.txt

# plot the results in a heatmap
python utils/hclust2.py -i output/merged_abundance_table.txt -o output_images/abundance_heatmap.png --skip_rows 1 --ftop 25 --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio .5 -l
--flabel_size 6 --legend_file HMP.log_scale.legend.png --max_flabel_len 200 --max_slabel_len 100  --slabel_size 6 --minv 0.1 --dpi 300
```
