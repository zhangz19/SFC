

rm(list=ls())
require(R.matlab)
source('util.R')

o <- readMat('datLM.mat')
Y <- o$Y
nTime <- nrow(Y)
nLoc <- ncol(Y)
nVar <- dim(o$X)[2]

#+++++++++++++++++++++++++++ kmeans with GAP stat: based on control run, no covariates
# http://www.sthda.com/english/wiki/print.php?id=239

require(cluster)
require(factoextra)
# 
# set.seed(20)
# ptm <- proc.time()[3]
# m0 <- clusGap(t(Y), kmeans, nstart = 25, K.max = 10, B = 500, verbose = TRUE)
# cputime <- as.numeric(proc.time()[3]-ptm)/60
# cat('\nCPUtime', round(cputime,3), 'minutes: completed!','\n') #CPUtime 15.93 minutes: completed! 
# save(file='m0_kmeansGap.Rdata', m0)
load('m0_kmeansGap.Rdata')
# # plot(m0, frame = FALSE, xlab = "Number of clusters k")

png('./tex/K_kmeans.png', height=3.4, width=6.4, units='in', pointsize=9, res=600)
fviz_gap_stat(m0, maxSE = list(method = "Tibs2001SEmax"))
dev.off()

print(m0, method = "Tibs2001SEmax")
# print(m0, method = "firstSEmax", SE.factor = 1)
K <- nrow(a <- m0[[1]])
kopt <- which(m0[[1]][1:(K-1),3] >= (m0[[1]][2:K,3] - m0[[1]][2:K,4]))[1] #Tibs2001SEmax
set.seed(11)
m1 <- kmeans(t(Y), centers=kopt, algorithm='Hartigan-Wong', nstart=25)
lab <- m1$cluster
lab1 <- sortLab(lab) # table(lab, lab1)
# save(file=paste0('lab_kmeans',kopt),lab1)
plotMapLab(lab1, Y, varids=1, saveFig=F, fnam=paste0('./tex/L_kmeans',kopt,'.png'))



kmax <- 8
set.seed(11)
m2 <- kmeans(t(Y), centers=kmax, algorithm='Hartigan-Wong', nstart=25)
lab <- m2$cluster
lab1 <- sortLab(lab)

source('util.R')
plotMapLab(lab1, Y, varids=1, saveFig=F, fnam=paste0('./tex/L_kmeans',kmax,'.png'))



W <- o$W
writeMat(paste0('lab_kmeans',kopt,'.mat'), lab1=as.numeric(lab1), W=W)




#+++++++++++++++++++++++++++ kmeans with GAP stat: based on control run, with 1 covariates
Y1 <- t(Y)
set.seed(20)
ptm <- proc.time()[3]
m0 <- clusGap(t(Y), kmeans, nstart = 25, K.max = 10, B = 500, verbose = TRUE)
cputime <- as.numeric(proc.time()[3]-ptm)/60
cat('\nCPUtime', round(cputime,3), 'minutes: completed!','\n') #CPUtime 15.93 minutes: completed!
save(file='m0_kmeansGap.Rdata', m0)





#+++++++++++++++++++++++++++ kmlcov
# install.packages('kmlcov')
require('kmlcov')
source('myglmClust.R') # modified the output format from their function

# # replicate their examples
# data(artifdata)
# write.csv(file='artifdata.csv',artifdata,row.names=F)

# res <- glmClust(formula = Y ~ clust(time + time2 + time3) + pop(treatTime), data = artifdata, ident = 'id', timeVar = 'time', effectVar = 'treatment', nClust = 4)
# 
# res <- kmlCov(formula = Y ~ clust(time + time2), data = artifdata, ident = 'id', timeVar = 'time', effectVar = 'treatment', nClust = 2:4, nRedraw = 2) #run 2 times for each cluster
# 
# a <- which_best(res, crit = "log-class-likelihood")  # best among nRedraw

# res <- myglmClust(formula = Y ~ clust(time + time2 + time3) + pop(treatTime), data = artifdata, ident = 'id', timeVar = 'time', effectVar = 'treatment', nClust = 4)
# 
# res$partition

# now try our data
dat <- data.frame(Y=as.vector(o$Y))
dat$time <- rep(1:nTime, nLoc)
dat$site <- rep(1:nLoc, each=nTime)
tmp <- numeric(); for(j in 1:nVar) tmp <- cbind(tmp, as.vector(o$X[,j,]))
dat <- cbind(dat, tmp)
names(dat) <-c('Y','time','site',paste0('X',1:nVar))
dat$time2 <- dat$time^2
dat$time3 <- dat$time^3

# run jobs with runkml.R

# set.seed(1)
# ptm <- proc.time()[3]
# m <- myglmClust(formula = Y ~ clust(time+time2+time3+X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15), data = dat, ident = 'site', timeVar = 'time', nClust = 6, max_itr=5000)
# cputime <- as.numeric(proc.time()[3]-ptm)/60
# cat('\nCPUtime', round(cputime,3), 'minutes: completed!','\n')  # CPUtime 7.619 minutes: completed!
# save(file='m_kmlcov.Rdata', m)
# load('m_kmlcov.Rdata')

# lab <- sortLab(m$partition)
# plotMapLab(lab, Y, varids=1, saveFig=FALSE, fnam=NULL)

load('sout_kml')
lab1 <- sortLab(lab) # table(lab, lab1)

png('./tex/K_kml.png', height=3.4, width=6.4, units='in', pointsize=12, res=600)
par(mar=c(2,2,1.2,.2)+0,mgp=c(1,.2,0), tck=-0.01, cex.axis=.85, cex.lab=1, cex.main=1)
plot(BIC~Ks, type='n', data=mat, xlim=c(2,11),ylab='BIC',xlab='Number of clusters k', main='Optimal number of clusters')
lines(BIC~Ks, data=mat,col='steelblue')
points(BIC~Ks, data=mat,col='steelblue',pch=16)
abline(v=kopt, lty=2, col='steelblue')
dev.off()

source('util.R')
plotMapLab(lab1, Y, varids=1, saveFig=T, fnam=paste0('./tex/L_kml',kopt,'.png'))



# not run




