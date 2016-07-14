# File: megan_csv_import.R
# Auth: u.niazi@imperial.ac.uk
# Date: 07/07/2016
# Desc: import data from megan csv output


if(!require(LearnBayes)) stop('R Package LearnBayes required')
## internal functions
# get alpha values for dirichlet posterior, using jeffery's non-informative prior
getAlpha = function(df, prior=1/2){
  seq = df$V2
  names(seq) = df$V1
  # remove the root node
  i = which(names(seq) == 'root')
  seq = seq[-i]
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
## data import
dfData = read.csv(file.choose(), stringsAsFactors = F, header=F)

ivAlpha = getAlpha(dfData)
# removed unassigned
i = which(names(ivAlpha) == 'Not assigned')
ivAlpha = ivAlpha[-i]
# simulate dirichlet posterior sample
mDir.post = getPosterior(ivAlpha)

## get the average vector and plot
iAve = colMeans(mDir.post)
head(iAve)

# break into groups
groups = cut(iAve, breaks = quantile(iAve, 0:20/20), include.lowest = T, labels = 1:20)

plot(sort(iAve[groups == '20']), type='l', main='Top 5% Abundant Taxa',
     xlab='Taxa', ylab='Proportion')

cvTop = names(sort(iAve[groups == '20'], decreasing = T))

## plot the top 20 samples
mPlot = mDir.post[,cvTop[1:20]]
plot.bar(mPlot, title = 'Sample 004_S1')

pie(colMeans(mPlot), cex=0.5, radius=1, angle=45)
