# File: metaphlan2_import.R
# Auth: u.niazi@imperial.ac.uk
# Date: 30/06/2016
# Desc: import data from metaphlan2 merged table


######################################
## data import

f_import.metaphlan = function(xls.file){
  # read the data in merged_table.txt file
  # it is a tab separated file with n lines of comments
  # the comment lines start with a #, so skip these lines
  infile = file(xls.file, 'rt')
  input = readLines(infile, n = 1000)
  n = grep('^#', x = input)
  close(infile)
  # we have n lines of comments, so skip n+1 lines  
  df = read.csv(file=xls.file, header=T, sep='\t', skip=length(n)+1,
                row.names=1)
  return(df)
}

dfData = f_import.metaphlan(file.choose())

## get genus level indices
rn = rownames(dfData)
i = grep('\\|g__\\w+$', rn, perl = T)

dfGenus = dfData[i,]
rn = rownames(dfGenus)
rn = strsplit(rn, '\\|')
rn = sapply(rn, function(x) x[length(x)])
rownames(dfGenus) = rn






strsplit(rn[2], '\\|')
grep('\\|g__\\w+$', head(rn), perl = T)
