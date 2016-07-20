# File: megan_multi_plots.R
# Auth: u.niazi@imperial.ac.uk
# Date: 20/07/2016
# Desc: import files from megan species output 

args = commandArgs(trailingOnly = TRUE)
csHelp = 'Run script with: \nRscript megan_multi_plots.R infile1.txt [infile2.txt ...]'

# check argument size
if (args[1] == 'h' || args[1] == 'help' || length(args) < 1) stop(cat(csHelp))

if(!require(LearnBayes)) stop('R Package LearnBayes required')
## internal functions
# get alpha values for dirichlet posterior, using jeffery's non-informative prior
getAlpha = function(df, prior=1/2){
  seq = df$V2
  names(seq) = df$V1
  # remove the root node
  i = which(names(seq) == 'root')
  # check if i is empty
  if (length(i) != 0) {seq = seq[-i]}
  alpha = seq + prior
  return(alpha)
}

# get posterior theta from posterior dirichlet
getPosterior = function(alpha, n=1000){
  p = rdirichlet(n, alpha)
  colnames(p) = names(alpha)
  #m = colMeans(p)
  return(p)
}

# bar plot with error bars
plot.bar = function(mDat, title='Abundance'){
  # get the median to plot
  p.old = par(mar=c(6,3,2,2)+0.1)
  mBar = apply(mDat, 2, mean)
  names(mBar) = colnames(mDat)
  yl = max(apply(mDat, 2, quantile, 0.98))
  l = barplot(mBar, beside=T, xaxt='n', ylim=c(0, yl), main=title)
  axis(side = 1, l[,1], labels=F)
  text(l[,1], y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=names(mBar), srt=45, adj=1, xpd=TRUE, cex=0.6)
  ## draw error bars
  f_barplot_errorbars = function(x.loc, y.loc){
    segments(x.loc, y.loc[1], x.loc, y.loc[2])
    segments(x.loc-0.1, y.loc[1], x.loc+0.1, y.loc[1])
    segments(x.loc-0.1, y.loc[2], x.loc+0.1, y.loc[2])
  }
  sapply(seq_along(1:ncol(mDat)), function(x) f_barplot_errorbars(l[x,1], quantile(mDat[,x], c(0.025, 0.975))))
  par(p.old)
}

######################################
## process one file at a time
sapply(args, function(x){
  # check if file exists
  if (!file.exists(x)) {
    stop(cat(paste('file', x, 'not found')))
  }
  # load the file
  dfData = read.csv(x, stringsAsFactors = F, header=F)
  # create short name for the file 
  cvName = unlist(strsplit(x, '/'))
  cvName = cvName[length(cvName)]
  # get posterior
  ivAlpha = getAlpha(dfData)
  # simulate dirichlet posterior sample
  mDir.post = getPosterior(ivAlpha)
  
  ## get the average vector and plot
  iAve = colMeans(mDir.post)
  
  # break into groups
  groups = cut(iAve, breaks = quantile(iAve, 0:10/10), include.lowest = T, labels = 1:10)
  
  # open device for plotting
  png(paste(x, '1.png', sep=''), width=609, height = 379, res=100)
  plot(sort(iAve[groups == '10']), type='l', main='Proportions for Top 10% Abundant Taxa',
       xlab='Taxa', ylab='Proportion', sub=cvName)
  dev.off(dev.cur())
  # get the top 10% sorted taxa
  cvTop = names(sort(iAve[groups == '10'], decreasing = T))
  
  ## plot the top 10% taxa
  mPlot = mDir.post[,cvTop]
  png(paste(x, '2.png', sep=''), width=609, height = 379, res=100)
  plot.bar(mPlot, title = paste('Top 10%', cvName))
  dev.off(dev.cur())
  
  png(paste(x, '3.png', sep=''), width=609, height = 379, res=100)
  pie(colMeans(mPlot), cex=0.7, radius=1, angle=45)
  dev.off(dev.cur())
  ## plot top 30% the taxa
  cvTop = names(sort(iAve[groups %in% c('8', '9', '10')], decreasing = T))
  
  ## plot the top 30% taxa
  mPlot = mDir.post[,cvTop]
  png(paste(x, '4.png', sep=''), width=609, height = 379, res=100)
  plot.bar(mPlot, title = paste('Top 30%', cvName))
  dev.off(dev.cur())
  png(paste(x, '5.png', sep=''), width=609, height = 379, res=100)
  pie(colMeans(mPlot), cex=0.5, radius=1, angle=45)
  dev.off(dev.cur())
})



