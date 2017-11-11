
# Codes atempt to inplementing ScFCM in Carlo et al. (2017)
# Here we consider balanced data with m_i = nt
require(splines)
require(nnet)

plotIt <- TRUE  # visual check during model fit

#================================================ prepare the data
o <- readMat('datLM.mat')
Y <- o$Y
X <- o$X[,1,]
W0 <- o$W0
nt <- nrow(Y)
ns <- ncol(Y)
time <- 1:nt

load('lab_kmeans7')  # find a good initial partition
d <- max(lab <- lab1)
# indr <- as.list(rep(NA), d); for(r in 1:d) indr[[r]] <- which(lab==r)
labmat <- matrix(lab[W0], nrow=ns)  # label matrix of neighbors

#+++++++++++++++ Choose K. the number of basis functions
getB <- function(y, k) return(bs(time, knots=seq(1,nt,len=k)[2:(k-1)], Boundary.knots=c(1,nt), degree=1, intercept=FALSE)) #degree=3 for cubic splines;  knots should be within the range

getBIC <- function(y, k){
  x <- getB(y, k)
  m <- lm(y ~ x ) 
  # plot(y~time); lines(fitted(m), col='red')
  return( nt*log(mean(m$residuals^2)) + k*log(nt) )
}

getBICk <- function(k) return(sum(apply(Y,2,getBIC, k=k)))

ks <- 1:45
bics <- sapply(ks, getBICk)
plot(bics ~ ks, type='b',pch=16)
kmin <- 17 #which.min(bics) #29

# design matrix for thetaY
mat <- lapply(seq_len(ncol(Y)), function(i) Y[,i])  # covert Y to list with column vectors
mat <- lapply(mat, getB, k=kmin)
mat <- do.call(rbind, mat)  # this is now the NT by k-1 design matrix
all(mat[(23-1)*nt+c(1:nt), ] == getB(Y[,23],kmin))  # a useful check: must be TRUE
mat <- as.data.frame(mat) 
names(mat) <- paste0('B',2:(1+ncol(mat)))  #1 corrsponds to intercept
mat$y <- as.vector(Y)

# design matrix for thetaZ
# count number of neighbors with label r
vat <- matrix(0, ns, d); for(r in 1:d) vat[,r] <- rowSums(labmat==r, na.rm=T)
vat <- cbind(vat, t(X))
vat <- data.frame(vat); names(vat) <- c( paste0('v',1:d), paste0('x1-',1:nt))
vat$y <- lab
m <- multinom(y~., data=vat)
thetaZ <- coef(m)

#+++++++++++++++ M-step
thetaY <- matrix(NA, d, 1+kmin)
labt <- rep(lab, each=nt)
if(plotIt) par(mfrow=c(2,4), mar=c(1,2,1.2,.2)+0,mgp=c(1,.2,0), tck=-0.01, cex.axis=.85, cex.lab=1, cex.main=1)
for(r in 1:d){
  indrt <- which(labt==r)
  m <- lm(y~., data=mat, subset=indrt)
  if(plotIt){ 
    yhat <- fitted(m)
    ym <- rowMeans(matrix(yhat, nrow=nt))
    matplot(matrix(mat$y[indrt], nrow=nt), type='l', lty=1, col='gray70', ylab='radiance', xlab='')
    lines(ym, col='red', lwd=2) 
  }
  coefs <- coef(m)
  thetaY[r,] <- c(coefs, sum(m$residuals^2) / (length(m$residuals) - length(coefs)))
}



#+++++++++++++++ R-step


# not run

