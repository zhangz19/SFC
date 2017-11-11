

repl <- as.numeric(commandArgs(TRUE))

require(R.matlab)
require(kmlcov)
source('myglmClust.R') # modified the output format from their function



o <- readMat('datLM.mat')
Y <- o$Y
nTime <- nrow(Y)
nLoc <- ncol(Y)
nVar <- dim(o$X)[2]

# now try our data
dat <- data.frame(Y=as.vector(o$Y))
dat$time <- rep(1:nTime, nLoc)
dat$site <- rep(1:nLoc, each=nTime)
tmp <- numeric(); for(j in 1:nVar) tmp <- cbind(tmp, as.vector(o$X[,j,]))
dat <- cbind(dat, tmp)
names(dat) <-c('Y','time','site',paste0('X',1:nVar))
dat$time2 <- dat$time^2
dat$time3 <- dat$time^3

K <- repl

set.seed(repl*10)
ptm <- proc.time()[3]
m <- myglmClust(formula = Y ~ clust(time+time2+time3+X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15), data = dat, ident = 'site', timeVar = 'time', nClust = K, max_itr=5000)
cputime <- as.numeric(proc.time()[3]-ptm)/60
cat('\nCPUtime', round(cputime,3), 'minutes: completed!','\n')

lab <- m$partition
cri <- m$criteria

save(file=paste0('out',repl), lab, cri, cputime)


combinIt <- FALSE
if(combinIt){
  Ks <- 2:20; cris <- numeric()
  for(repl in Ks){
    load(paste0('out',repl))
    cris <- rbind(cris, cri)
  }
  mat <- data.frame(cbind(Ks, cris))
  y <- diff(mat$BIC)
  kopt <- mat$Ks[which(y>0)[1]]  #optimal K
  load(paste0('out',kopt))
  save(file='sout_kml', mat, lab, kopt)  #lab=label under optimal partition
}





# not run


