

#++++++++++++++++ 2017-10-10 prepare datLM.mat and some graphics
rm(list=ls())

require(maptools)
library(maps)
library(fields)
library(lattice)
library(ggplot2)
#library(ape)
library(R.matlab)
source('util.R')

dirs <- './LatMonData/'
load(paste(dirs,'LatMon.RD',sep='')) # from other methods

# rad2.diff.mon[wnumber,lat,covariates,month]
#wave number=50:2760 cm^-1
#lat = seq(-85,85 by=10) # -85S:85N
#covariates: 1- control run  2:16- 15 different kernels
#month=200301-200412
# 2014-5-13: use original month
nx <- length(month <- c(1:24)) #/24
# nx <- length(month <- c(0:23)/23)
# 2013-11-5 also map lattitude into (0,1)
# ny <- length(lat <- seq(-85,85, by=10)/360*2*pi)
lat <- seq(-85,85, by=10)#/360*2*pi;  # 2014-5-13: use degree
# 2014-5-13: no need to transform to [0,1]
# lat <- (lat-min(lat))/(max(lat)-min(lat))  # transform to [0,1]
ny <- length(lat)

# tnf <- length(induse <- 1:dim(rad2.diff.mon)[1] )
# # (tnf <- length(induse <- c(153:2200)))
# 
# matplot(rad2.diff.mon[induse,,1,1], pch=16, type='b', cex=.6)

#matplot(rad2.diff.mon[1:1700,,1,1], pch=16, type='b', cex=.6)


useAllPoints <- FALSE #True use 
if(!useAllPoints){
  nw <- 1700 #dim(rad2.diff.mon)[1] #1700 #dim(rad2.diff.mon)[1]           #2e3; #
  tnf <- 128  #log2(nw)
  induse <- seq(1, nw, by=round(nw/tnf)); induse <- induse[1:tnf]
}
if(useAllPoints) induse <- 1:dim(rad2.diff.mon)[1]  # to use the whole curve
##inds <- round(seq(1, nw, length=tnf)); inds <- inds[1:tnf]
y <- rad2.diff.mon[induse,,1,1]
colMeans(y)[1:5]
# matplot(y, type='l', lty=1,pch=16)

# check with table 1 in Kato 2011
# a <- apply(rad2.diff.mon, 2, range)


#vec <- 1:1700
#tmpind <- seq(1, nw, by=round(1700/2^7));
# print(round(1700/2^7))


# 2015-5-11
# jpeg('./tex/argonneFig1.jpg',quality=100, height=4000,width=6000, pointsize=16, res=600)

# mat <- rad2.diff.mon[inds,,1,]

# reshape data to T by p by N
plotit <- F
dat <- array(NA, dim=c(tnf, dim(rad2.diff.mon)[3], nx*ny))
for (varid in 1:dim(rad2.diff.mon)[3]) {    #dim(rad2.diff.mon)[3]
  ran <- range(rad2.diff.mon[induse, , varid, ])
  scal1 <- 80; scal2 <- .0005; 
  scal3 <- max(rad2.diff.mon[induse, , 1, ])
  k <- 1
  if(plotit){
    # matplot(dat,type='b',pch=16,col='gray50')
    # palette(rainbow(d))
    par(mfrow=c(1,1), mar=c(2,2.5,0,0)+.3,mgp=c(1.1,.2,0), tck=-0.01, cex.axis=.8, cex.lab=1.2, cex.main=1, xaxt='n', yaxt='n')
    plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude',main=''); #grid(col='gray70')     main=paste('variable id =',varid) Spectral radiance difference over 24 months (2003/01-2004/12)
  }
  for(i in 1:nx){ #month
    for(j in 1:ny){ #latitude first
      dat[,varid,k] <- y <- rad2.diff.mon[induse, j, varid, i]
      y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
      x1 <- c(1:length(y1))*scal2+month[i]
      if(plotit) lines(y1 ~ x1, lwd=1, col='gray20') #,d+1-labmat[k,varid]
      #bp <- round(Box.test(y, type="Ljung-Box", lag=1)$p.value,3)
      #text(x1[1],y1[1]-.01,labels=as.character(bp),col='gray20',cex=.8)
      k <- k+1
    }
  }
  # write.table(dat, file=paste('dat',varid,'.txt',sep=''),row.names=F,col.names=F)
  if(plotit){
    par(xaxt='s', yaxt='s')
    axis(side=1,at=1:24,labels=seq(1,24,by=1))
    par(las=1,mgp=c(1.7,.5,0), tck=-0.01); axis(side=2,at=lat,labels=lat)
  }
}


# dev.off()


# extract neighrbor matrix
nmonth <- 1:nx
nlat <- 1:ny
n <- nx*ny
coords <- matrix(0, n, 2) # only used for calculating the Manhattan distance below
W <- matrix(0, n, n)
k <- 1
for(i in 1:nx){
  for(j in 1:ny){
    coords[k,] <- c(nmonth[i], nlat[j])
    k <- k+1
  }
}


# # just use Manhattan  --- this is for 4-NN
# W <- as.matrix(dist(coords, method='manhattan'))
# range(rowSums(W==1))
# W0 <- matrix(n+1, n, max(rowSums(W==1)))  # append "n+1" in the end
# for(i in 1:n){
#   ind <- which(W[i,]==1)
#   W0[i, 1:length(ind)] <- ind
# }

# use 8-NN which seems more proper in this case with irregular shape of clusters
W <- ( as.matrix(dist(coords, method='euclidean')) <= sqrt(2)+1e-10 )*1
W <- W + diag(nrow(W))
W0 <- matrix(n+1, n, max(rowSums(W==1)))  # append "n+1" in the end
for(i in 1:n){
  ind <- which(W[i,]==1)
  W0[i, 1:length(ind)] <- ind
}


# now calculate the coords in raw variable values to save
coords <- matrix(0, n, 2)
k <- 1
for(i in 1:nx){
  for(j in 1:ny){
    coords[k,] <- c(month[i], lat[j])
    k <- k+1
  }
}


# 2017-10-20: now put all data processing here, and the data will not change any more
# dat[,-1,] <- 1e4*dat[,-1,]  # scale transform
# 10-27: how about scale by the range? 
for(i in 1:dim(dat)[2]) dat[,i,] <- dat[,i,]/(max(dat[,i,]) - min(dat[,i,]))

Y <- dat[,1,]  #response for control run
# matplot(Y, type='l', lty=1)
X <- dat[,-1,]  # these are functional covariates, the radience change under perturbed runs



# save results: D is the distance matrix, W is the adjacency matrix, W0 is NbyK storing Kth nearest neighbors, Y response, X covariates, coords: Nby2 for 2-d input space
writeMat('datLM.mat', D=W, W=(W==1)*1.0, W0=W0, Y=Y, X=X, coords=coords)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ get figures for the paper
Ynams <- c('Control run', 'Near-surface&skin temperature', 'Surface-200-hPa temperature', '200-10-hPa temperature', 'Surface-500-hPa water vapor', '500-200-hPa water vapor', 'Low-level cloud-top height', 'Midlevel cloud-top height', 'High-level cloud-top height', 'Low-level cloud fraction', 'Midlevel cloud fraction', 'High-level cloud fraction', 'Thin ice cloud optical thickness', 'Ice cloud particle size', 'Water cloud optical thickness', 'Water cloud particle size')


pchs <- c(17,15,16,2,22,25)
cols0 <- c('#4575b4','#74add1','#abd9e9','#e0f3f8','#ffffbf','#fee090','#fdae61','#f46d43','#d73027','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788') #latitude
cols0 <- rev(cols0) #dark purple: south pole; dark blue: north pole
cols <- rep(cols0, nx)
#plot(1:ny, col=cols0)



png(paste0('./tex/All.png'), height=5,width=7.1, units='in', pointsize=12, res=600)

m0 <- matrix(c(1:16,rep(17,4)),nrow = 5,ncol = 4,byrow = TRUE)
layout(mat = m0, heights = c(rep(.22,4),1-(.22*4)))
par(oma=c(1.2,1.2,.1,0), font.lab=1, cex.lab=1.1,cex.axis=.9, cex.main=.95, mar=c(1,1,1,0)+.5, mgp=c(1,.3,0),tck=-.02,las=1,xaxt='n')

# par(mfrow=c(4,4), mar=c(2,2,1,.3)+.5,mgp=c(1.2,.3,0), tck=-0.01, cex.axis=1, cex.lab=1, cex.main=1)
nt <- length(induse)
vec <- seq(1,nt,by=round(nt/6))
for(j in 1:length(Ynams)){
  matplot(dat[,j,], type='l', lty=1, col=cols, main=Ynams[j], xlab='', ylab='')
  abline(h=0, col='gray50', lwd=1, lty=1)
  par(xaxt='s')
  axis(side=1,at=vec,labels=FALSE,las=1)
  text(x=vec, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), labels=induse[vec], srt=15, adj=1,xpd=T,cex=.9)
  par(xaxt='n',xpd=F) 
}


source('util.R')
latNam <- lat2lab(lat, padding=T)
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend('center',legend=latNam, bty='n', fill=NA, border=NA, text.col='black',inset=0, ncol=6, cex=1, pt.cex=20, pch='.',col=cols0)

par(las=0)
mtext(c(
  # as.expression(substitute(paste(a1,"",mu,a2),list(a1="Compound concentration (", a2="g/L)") )),
  paste("Wavenumber", parse(text="(cm^-1)")),
  "Spectral radiance change"
  ), side=1:2, line=c(.3,.2), outer=TRUE, cex=.8)

dev.off()









#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ backup chunk

# investigate latitude effect
mat <- as.vector(dat)
mat <- as.data.frame(mat); names(mat) <- 'y'
mat$lat <- rep(rep(lat, each=tnf), nx)
mat$lat <- as.factor(mat$lat)
mat$month <- rep(month, each=ny*tnf)
mat$x <- rep(1:tnf, nx*ny)
mat$id <- rep(1:(nx*ny), each=tnf)

# 2015-3-11: get labels of lat/mon for the 432 sites
lats <- rep(lat, nx)
mons <- rep(month, each=ny)

# writeMat('./inds.mat',lats=lats, mons=mons)

# a simple check
matplot(dat[,lats<=-70], pch=16,col=1)


lattice.options(default.theme = canonical.theme(color = FALSE))
my.padding <- list(strip = .8, top.padding = 0, main.key.padding = 0, 
                   key.axis.padding = 0,
                   axis.ylab.padding = 0, ylab.key.padding  = 0, right.padding = 0,
                   left.padding = 0, key.sub.padding = 0, bottom.padding = 0,
                   axis.right = 0, key.right = 0)

xyplot(y ~ x|lat, type='l', group=id, data=mat, xlab='Wave number', main='', scale=list(y=list(cex=.6, rot=c(45,45),tck=c(.2,.2)), x=list(cex=.6, rot=c(0,0),tck=c(.5,.5))), #relation='free', 
       strip = strip.custom(par.strip.text = list(cex=.8, font=3)),  
       # panel=function(x,y,...){panel.xyplot(x,y, pch=16, cex=.5);  
       # if(var(x, na.rm=T)!=0) panel.abline(lm(y~x), col='red')
       # }, 
       par.settings = list(layout.heights = my.padding, layout.widths = my.padding) #, layout=c(6,8)
)


xyplot(y ~ x|month, type='l', group=id, data=mat, xlab='Wave number', main='', scale=list(y=list(cex=.6, rot=c(45,45),tck=c(.2,.2)), x=list(cex=.6, rot=c(0,0),tck=c(.5,.5))), #relation='free', 
       strip = strip.custom(par.strip.text = list(cex=.8, font=3)),  
       # panel=function(x,y,...){panel.xyplot(x,y, pch=16, cex=.5);  
       # if(var(x, na.rm=T)!=0) panel.abline(lm(y~x), col='red')
       # }, 
       par.settings = list(layout.heights = my.padding, layout.widths = my.padding) #, layout=c(6,8)
)







# 2/13/2015: final data processing for complete data
rm(list=ls())


dirs <- './LatMonData/'
load(paste(dirs,'LatMon.RD',sep='')) # from other methods

# rad2.diff.mon, 2711   18   16   24



# rad2.diff.mon[wnumber,lat,covariates,month]
#wave number=50:2760 cm^-1
#lat = seq(-85,85 by=10) # -85S:85N
#covariates: 1- control run  2:16- 15 different kernels
#month=200301-200412
# 2014-5-13: use original month
nx <- length(month <- c(1:24)) #/24
# nx <- length(month <- c(0:23)/23)
# 2013-11-5 also map lattitude into (0,1)
# ny <- length(lat <- seq(-85,85, by=10)/360*2*pi)
lat <- seq(-85,85, by=10)#/360*2*pi;  # 2014-5-13: use degree
# 2014-5-13: no need to transform to [0,1]
# lat <- (lat-min(lat))/(max(lat)-min(lat))  # transform to [0,1]
ny <- length(lat)

tnf <- length(induse <- 1:dim(rad2.diff.mon)[1] )
#(tnf <- length(induse <- c(153:2200)))

matplot(rad2.diff.mon[induse,,1,1], pch=16, type='b', cex=.6)

#matplot(rad2.diff.mon[1:1700,,1,1], pch=16, type='b', cex=.6)

#nw <- dim(rad2.diff.mon)[1] #1700 #dim(rad2.diff.mon)[1]           #2e3; #
#tnf <- log2(tlen)
#induse <- seq(1, nw, by=round(nw/tnf)); induse <- induse[1:tnf]
##inds <- round(seq(1, nw, length=tnf)); inds <- inds[1:tnf]
#matplot(rad2.diff.mon[induse,,1,1], type='l', lty=1,pch=16)

#vec <- 1:1700
#tmpind <- seq(1, nw, by=round(1700/2^7));
print(round(1700/2^7))


# 2015-5-11
jpeg('argonneFig1.jpg',quality=100, height=4000,width=6000, pointsize=16, res=600)

# mat <- rad2.diff.mon[inds,,1,]

# reshape data
dat <- array(NA, dim=c(tnf, dim(rad2.diff.mon)[3], nx*ny))
for (varid in 1:dim(rad2.diff.mon)[3]) {    #dim(rad2.diff.mon)[3]
  ran <- range(rad2.diff.mon[induse, , varid, ])
  scal1 <- 80; scal2 <- .0005; 
  scal3 <- max(rad2.diff.mon[induse, , 1, ])
  plotit <- FALSE
  k <- 1
  if(plotit){
    # matplot(dat,type='b',pch=16,col='gray50')
    # palette(rainbow(d))
    par(mfrow=c(1,1), mar=c(2,2.5,0,0)+.3,mgp=c(1.1,.2,0), tck=-0.01, cex.axis=.8, cex.lab=1.2, cex.main=1, xaxt='n', yaxt='n')
    plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude',main=''); #grid(col='gray70')     main=paste('variable id =',varid) Spectral radiance difference over 24 months (2003/01-2004/12)
  }
  for(i in 1:nx){
    for(j in 1:ny){
      dat[,varid,k] <- y <- rad2.diff.mon[induse, j, varid, i]
      y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
      x1 <- c(1:length(y1))*scal2+month[i]
      if(plotit) lines(y1 ~ x1, lwd=1, col='gray20') #,d+1-labmat[k,varid]
      #bp <- round(Box.test(y, type="Ljung-Box", lag=1)$p.value,3)
      #text(x1[1],y1[1]-.01,labels=as.character(bp),col='gray20',cex=.8)
      k <- k+1
    }
  }
  # write.table(dat, file=paste('dat',varid,'.txt',sep=''),row.names=F,col.names=F)
  if(plotit){
    par(xaxt='s', yaxt='s')
    axis(side=1,at=1:24,labels=seq(1,24,by=1))
    par(las=1,mgp=c(1.7,.5,0), tck=-0.01); axis(side=2,at=lat,labels=lat)
  }
}


dev.off()


library(R.matlab)
writeMat('datComplete.mat',Y=dat[,1,], X=dat[,-1,])





# investigate latitude effect
mat <- as.vector(dat)
mat <- as.data.frame(mat); names(mat) <- 'y'
mat$lat <- rep(rep(lat, each=tnf), nx)
mat$lat <- as.factor(mat$lat)
mat$month <- rep(month, each=ny*tnf)
mat$x <- rep(1:tnf, nx*ny)
mat$id <- rep(1:(nx*ny), each=tnf)

# 2015-3-11: get labels of lat/mon for the 432 sites
lats <- rep(lat, nx)
mons <- rep(month, each=ny)

writeMat('../inds.mat',lats=lats, mons=mons)

# a simple check
matplot(dat[,lats<=-70], pch=16,col=1)


lattice.options(default.theme = canonical.theme(color = FALSE))
my.padding <- list(strip = .8, top.padding = 0, main.key.padding = 0, 
                   key.axis.padding = 0,
                   axis.ylab.padding = 0, ylab.key.padding  = 0, right.padding = 0,
                   left.padding = 0, key.sub.padding = 0, bottom.padding = 0,
                   axis.right = 0, key.right = 0)

xyplot(y ~ x|lat, type='l', group=id, data=mat, xlab='Wave number', main='', scale=list(y=list(cex=.6, rot=c(45,45),tck=c(.2,.2)), x=list(cex=.6, rot=c(0,0),tck=c(.5,.5))), #relation='free', 
       strip = strip.custom(par.strip.text = list(cex=.8, font=3)),  
       # panel=function(x,y,...){panel.xyplot(x,y, pch=16, cex=.5);  
       # if(var(x, na.rm=T)!=0) panel.abline(lm(y~x), col='red')
       # }, 
       par.settings = list(layout.heights = my.padding, layout.widths = my.padding) #, layout=c(6,8)
)


xyplot(y ~ x|month, type='l', group=id, data=mat, xlab='Wave number', main='', scale=list(y=list(cex=.6, rot=c(45,45),tck=c(.2,.2)), x=list(cex=.6, rot=c(0,0),tck=c(.5,.5))), #relation='free', 
       strip = strip.custom(par.strip.text = list(cex=.8, font=3)),  
       # panel=function(x,y,...){panel.xyplot(x,y, pch=16, cex=.5);  
       # if(var(x, na.rm=T)!=0) panel.abline(lm(y~x), col='red')
       # }, 
       par.settings = list(layout.heights = my.padding, layout.widths = my.padding) #, layout=c(6,8)
)



# extract neighrbor matrix
nmonth <- 1:nx
nlat <- 1:ny
n <- nx*ny
coords <- matrix(0, n, 2)
W <- matrix(0, n, n)
k <- 1
for(i in 1:nx){
  for(j in 1:ny){
    coords[k,] <- c(nmonth[i], nlat[j])
    k <- k+1
  }
}

# just use Manhattan
W <- as.matrix(dist(coords, method='manhattan'))
range(rowSums(W==1))
W0 <- matrix(0, n, max(rowSums(W==1)))
for(i in 1:n){
  ind <- which(W[i,]==1)
  W0[i, 1:length(ind)] <- ind
}



# writeMat('../../codes2/datLM.mat', D=W, W=(W==1)*1.0, W0=W0, y=dat)



# investigate variance versus mean
ra <- 5
N <- nrow(mat)
tmp <- c(rep(NA, ra), mat$y, rep(NA, ra) )
vars <- numeric(N)
for(i in 1:N){
  tmpi <- i + seq(-ra,ra,by=1)+ra; #tmpi <- tmpi[-which(tmpi==0)]
  vars[i] <- var(mat$y[tmpi], na.rm=T)
}

plot(vars ~ mat$y) # irregular



# investigate variance versus wavenumber
vars <- numeric(nrow(mat))
for(i in 1:tnf){
  inds <- which(mat$x == i)
  vars[inds] <- var(mat$y[inds])
}

plot(vars ~ mat$x, type='l', lty=1,pch=16)



errs <- numeric(nrow(mat))
for(j in lat){
  for(i in 1:tnf){
    inds <- which(mat$x == i & mat$lat == j)
    errs[inds] <- var(mat$y[inds])
  }
}

xyplot(errs ~ x|lat, type='l', group=id, data=mat, xlab='Wave number', main='', scale=list(y=list(cex=.6, rot=c(45,45),tck=c(.2,.2)), x=list(cex=.6, rot=c(0,0),tck=c(.5,.5))), #relation='free', 
       strip = strip.custom(par.strip.text = list(cex=.8, font=3)),  
       # panel=function(x,y,...){panel.xyplot(x,y, pch=16, cex=.5);  
       # if(var(x, na.rm=T)!=0) panel.abline(lm(y~x), col='red')
       # }, 
       par.settings = list(layout.heights = my.padding, layout.widths = my.padding) #, layout=c(6,8)
)









mat$vars <- vars
xyplot(vars ~ y|x, type='p', data=mat, xlab='Wave number', main='', scale=list(y=list(cex=.6, rot=c(45,45),tck=c(.2,.2)), x=list(cex=.6, rot=c(0,0),tck=c(.5,.5))), #relation='free', 
       strip = strip.custom(par.strip.text = list(cex=.8, font=3)),  
       # panel=function(x,y,...){panel.xyplot(x,y, pch=16, cex=.5);  
       # if(var(x, na.rm=T)!=0) panel.abline(lm(y~x), col='red')
       # }, 
       par.settings = list(layout.heights = my.padding, layout.widths = my.padding) #, layout=c(6,8)
)




# test covariates effects functions
# reshape data

#X <- matrix(NA, nrow(mat), dim(rad2.diff.mon)[3]-1)
X <- array(NA, dim=c(nrow(dat), dim(rad2.diff.mon)[3]-1, ncol(dat)))
for (varid in 2:dim(rad2.diff.mon)[3]) { 
  print(varid)
  #ran <- range(rad2.diff.mon[inds, , varid, ])
  dat1 <- matrix(NA, tnf, nx*ny)
  #scal1 <- 80; scal2 <- .006; 
  #scal3 <- max(rad2.diff.mon[inds, , 1, ])
  k <- 1
  # palette(rainbow(d))
  #par(mfrow=c(1,1), mar=c(2,2,1,0)+.3,mgp=c(1.1,.2,0), tck=-0.01, cex.axis=.8, cex.lab=.8, cex.main=1)
  #plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude',main=paste('variable id =',varid)); #grid(col='gray70')
  for(i in 1:nx){
    for(j in 1:ny){
      dat1[,k] <- y <- rad2.diff.mon[induse, j, varid, i]
      # y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
      # x1 <- c(1:length(y1))*scal2+month[i]
      # lines(y1 ~ x1, lwd=1, col='gray20') #,d+1-labmat[k,varid]
      #bp <- round(Box.test(y, type="Ljung-Box", lag=1)$p.value,3)
      #text(x1[1],y1[1]-.01,labels=as.character(bp),col='gray20',cex=.8)
      k <- k+1
    }
  }
  # write.table(dat, file=paste('dat',varid,'.txt',sep=''),row.names=F,col.names=F)
  # X[,varid-1] <- as.vector(dat1)
  X[,varid-1,] <- dat1
}




coords <- matrix(0, n, 2)
k <- 1
for(i in 1:nx){
  for(j in 1:ny){
    coords[k,] <- c(month[i], lat[j])
    k <- k+1
  }
}
writeMat('../coords.mat', coords=coords)


writeMat('../datLM.mat', D=W, W=(W==1)*1.0, W0=W0, Y=dat, X=X)



#regression with all covariates, pooled at each wavenumber
tmp <- mat$y
errs <- numeric(nrow(mat))
for(i in 1:tnf){
  inds <- which(mat$x == i)
  errs[inds] <- resid(lm(mat$y[inds] ~ X[inds,]))
}

matplot(matrix(errs, nrow=tnf), type='l', lty=1,pch=16)





#regression with all covariates, pooled at each wavenumber, by lat
tmp <- mat$y
errs <- numeric(nrow(mat))
for(j in lat){
  for(i in 1:tnf){
    inds <- which(mat$x == i & mat$lat == j)
    # errs[inds] <- resid(lm(mat$y[inds] ~ X[inds,]))
    errs[inds] <- coef(lm(mat$y[inds] ~ X[inds,]))[2]
  }
}

xyplot(errs ~ x|lat, type='l', group=id, data=mat, xlab='Wave number', main='', scale=list(y=list(cex=.6, rot=c(45,45),tck=c(.2,.2)), x=list(cex=.6, rot=c(0,0),tck=c(.5,.5))), #relation='free', 
       strip = strip.custom(par.strip.text = list(cex=.8, font=3)),  
       # panel=function(x,y,...){panel.xyplot(x,y, pch=16, cex=.5);  
       # if(var(x, na.rm=T)!=0) panel.abline(lm(y~x), col='red')
       # }, 
       par.settings = list(layout.heights = my.padding, layout.widths = my.padding) #, layout=c(6,8)
)





#regression with all covariates, pooled at each wavenumber, by month
errs <- numeric(nrow(mat))
for(j in month){
  for(i in 1:tnf){
    inds <- which(mat$x == i & mat$month == j)
    # errs[inds] <- resid(lm(mat$y[inds] ~ X[inds,]))
    errs[inds] <- coef(lm(mat$y[inds] ~ X[inds,]))[2]
  }
}

xyplot(errs ~ x|month, type='l', group=id, data=mat, xlab='Wave number', main='', scale=list(y=list(cex=.6, rot=c(45,45),tck=c(.2,.2)), x=list(cex=.6, rot=c(0,0),tck=c(.5,.5))), #relation='free', 
       strip = strip.custom(par.strip.text = list(cex=.8, font=3)),  
       # panel=function(x,y,...){panel.xyplot(x,y, pch=16, cex=.5);  
       # if(var(x, na.rm=T)!=0) panel.abline(lm(y~x), col='red')
       # }, 
       par.settings = list(layout.heights = my.padding, layout.widths = my.padding) #, layout=c(6,8)
)




















print('done')

#for(i in 1:(n-1)){
#  for(j in (i+1):n){ 
#    W[i,j] <- sqrt(sum((coords[i,]-coords[j,])^2))    
#    W[j,i] <- W[i,j]
#  }
#}
#diag(W) <- 10
#W0 <- (W<=1)
#rowSums(W0)





#write.table(W, file='neighbor.txt',row.names=F,col.names=F)

#coords <- matrix(0, n, 2)
#k <- 1
#for(i in 1:nx){
# for(j in 1:ny){
#    coords[k,] <- c(month[i], lat[j])
#    k <- k+1
# }
#}
#write.table(coords, file='coords.txt',row.names=F,col.names=F)

# save distance vector
mydist <- numeric(n*(n-1)/2); k <- 1
for(i in 1:(n-1)){
  for(j in (i+1):n){ 
    mydist[k] <- sqrt(sum((coords[i,]-coords[j,])^2)); k <- k+1
  }
}
mydist <- round(mydist, 8)
L <- length(hvec <- sort(unique(mydist)))
write.table(mydist, file='mydist.txt',row.names=F,col.names=F)

# save distance matrix
mydist <- matrix(0,ncol=n, nrow=n)
for(i in 1:(n-1)){
  for(j in (i+1):n){ 
    mydist[i,j] <- sqrt(sum((coords[i,]-coords[j,])^2))
    mydist[j,i] <- mydist[i,j]
  }
}
write.table(mydist, file='mydist2.txt',row.names=F,col.names=F)


## 2014-5-13  save datLatMon2
a <- readMat('../../codes/datLatMon.mat')
Y <- a$Y
writeMat('datLatMon2.mat', coords=as.matrix(coords), Y=Y)













## backup

#----------------------------------- spatio-temporal model
rm(list=ls())
load('LatMon.RD') # from other methods
labmat <- readMat('../../codes/labmatLM.mat')$labmat
d <- 8 
# rad2.diff.mon[wnumber,lat,covariates,month]
#wave number=50:2760 cm^-1
#lat = seq(-85,85 by=10) # -85S:85N
#covariates: 1- control run  2:16- 15 different kernels
#month=200301-200412
# 2014-5-13: use original month
nx <- length(month <- c(1:24)) #/24
# nx <- length(month <- c(0:23)/23)
# 2013-11-5 also map lattitude into (0,1)
# ny <- length(lat <- seq(-85,85, by=10)/360*2*pi)
lat <- seq(-85,85, by=10)#/360*2*pi;  # 2014-5-13: use degree
# 2014-5-13: no need to transform to [0,1]
# lat <- (lat-min(lat))/(max(lat)-min(lat))  # transform to [0,1]
ny <- length(lat)

matplot(rad2.diff.mon[1:1700,,1,1], pch=16)

nw <- 1700 #dim(rad2.diff.mon)[1]           #2e3; #
tnf <- 2^7
inds <- seq(1, nw, by=round(nw/tnf)); inds <- inds[1:tnf]
#inds <- round(seq(1, nw, length=tnf)); inds <- inds[1:tnf]
matplot(rad2.diff.mon[inds,,1,1], type='l', lty=1,pch=16)

# reshape data
pdf('ProfilePlot.pdf', paper='a4r', width = 12, height = 8)
for (varid in 1:dim(rad2.diff.mon)[3]) {
  ran <- range(rad2.diff.mon[inds, , varid, ])
  dat <- matrix(NA, tnf, nx*ny)
  scal1 <- 1; scal2 <- .0003; scal3 <- max(rad2.diff.mon[inds, , 1, ])
  k <- 1
  palette(rainbow(d))
  plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude',main=paste('variable id =',varid)); #grid(col='gray70')
  for(i in 1:nx){
    for(j in 1:ny){
      dat[,k] <- y <- rad2.diff.mon[inds, j, varid, i]
      y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
      x1 <- c(1:length(y1))*scal2+month[i]
      # lines(y1 ~ x1, lwd=1, col=d+1-labmat[k,varid]) #col='gray20',
      #bp <- round(Box.test(y, type="Ljung-Box", lag=1)$p.value,3)
      #text(x1[1],y1[1]-.01,labels=as.character(bp),col='gray20',cex=.8)
      k <- k+1
    }
  }
  write.table(dat, file=paste('dat',varid,'.txt',sep=''),row.names=F,col.names=F)
}
dev.off()
# matplot(dat,type='b',pch=16,col='gray50')

nmonth <- 1:nx
nlat <- 1:ny
n <- nx*ny
coords <- matrix(0, n, 2)
W <- matrix(0, n, n)
k <- 1
for(i in 1:nx){
  for(j in 1:ny){
    coords[k,] <- c(nmonth[i], nlat[j])
    k <- k+1
  }
}
for(i in 1:(n-1)){
  for(j in (i+1):n){ 
    W[i,j] <- sqrt(sum((coords[i,]-coords[j,])^2))    
    W[j,i] <- W[i,j]
  }
}
write.table(W, file='neighbor.txt',row.names=F,col.names=F)

coords <- matrix(0, n, 2)
k <- 1
for(i in 1:nx){
  for(j in 1:ny){
    coords[k,] <- c(month[i], lat[j])
    k <- k+1
  }
}
write.table(coords, file='coords.txt',row.names=F,col.names=F)

# save distance vector
mydist <- numeric(n*(n-1)/2); k <- 1
for(i in 1:(n-1)){
  for(j in (i+1):n){ 
    mydist[k] <- sqrt(sum((coords[i,]-coords[j,])^2)); k <- k+1
  }
}
mydist <- round(mydist, 8)
L <- length(hvec <- sort(unique(mydist)))
write.table(mydist, file='mydist.txt',row.names=F,col.names=F)

# save distance matrix
mydist <- matrix(0,ncol=n, nrow=n)
for(i in 1:(n-1)){
  for(j in (i+1):n){ 
    mydist[i,j] <- sqrt(sum((coords[i,]-coords[j,])^2))
    mydist[j,i] <- mydist[i,j]
  }
}
write.table(mydist, file='mydist2.txt',row.names=F,col.names=F)


## 2014-5-13  save datLatMon2
a <- readMat('../../codes/datLatMon.mat')
Y <- a$Y
writeMat('datLatMon2.mat', coords=as.matrix(coords), Y=Y)





#------------------------- investigate whether we should use random effects
a <- readMat('../../codes/datLatMon.mat')
Y <- a$Y
N <- ncol(Y); T = nrow(Y)
par(mar=c(2,2,0,0)); matplot(Y, type='n'); for(s in 1:ncol(Y)) lines(1:T, Y[,s], col='gray80')

par(mar=c(2,2,0,0)); plot(1:T, vars <- apply(Y,1,var), pch=16, type='b'); lines(1:T, predict(loess(vars~c(1:T))), col='red', lwd=2);# lines(fit.sp <- spline(vars~c(1:T),n=T/2), col='blue', lwd=2)

hist(vars)















# 1-15-2014
varid <- 1
cols <- rep(1:18, 24)
ran <- range(Y)
scal1 <- .44; scal2 <- .00025; scal3 <- max(Y)
k <- 1
palette(rainbow(18))
par(mar=c(2.4,2.4,0,0)+.05,mgp = c(1.3, .2, 0),tck=-0.005,cex.axis=1,cex.lab=1, xaxt='n', yaxt='n')
plot(1, type='n', ylim=c(min(lat), max(lat)+.01), xlim=range(month),xlab='month',ylab='latitude',main=''); #grid(col='gray70')
par(xaxt='s', yaxt='s'); axis(1,at=month, label=1:24)
axis(2,at=lat, label=seq(-85,85,len=length(lat)))
for(i in 1:nx){
  for(j in 1:ny){
    y <- Y[,k]
    y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
    x1 <- c(1:length(y1))*scal2+month[i]
    lines(y1 ~ x1, lwd=1, col=cols[k]) #col='gray20',
    k <- k+1
  }
}













# for simulation, pick centroids
cts <- as.integer(round(quantile(1:N,seq(0,1,len=15))))
length(cts <- cts[-c(1,length(cts))])
labs <- numeric(N); for(i in 1:N) labs[i] <- which(mydist[i,cts] == min(mydist[i,cts]))[1]
d <- max(labs)
simud <- max(labs)
# plot(coords, type='n'); text(coords, labels=c(1:N))
plot(coords, type='n'); palette(rainbow(length(unique(labs)))); text(coords, labels=labs, col=labs); points(coords[cts,], pch=16, col='red', cex=1.2)
simuMu <- matrix(0,T,d)
for(r in 1:d) simuMu[,r] <- rowMeans(Y[,labs==r])
# simuMu <- simuMu[,-c(1:2)]
matplot(simuMu, pch=16, type='l', lwd=2, lty=1)
writeMat('../../codes/simusetup.mat',simuMu=simuMu, labs=labs, center=cts, d=d)

write.table(simuMu, file='simuMu.txt',row.names=F,col.names=F)


simus <- readMat('../../codes/FSCsimu7.mat')
waveT <- readMat('../../codes/X.mat')$X
Y0 <- simus$Y
labs <- as.numeric(simus$truePara[[3]]); d <- max(labs)




# -----------------------------------------------------------------------------------

# get Figure rdo
postscript(file='realdataonemonth.eps',pointsize=25,width=16,height=8,horizontal=F,paper='special')
lat0 <- seq(-85,85, by=10)
par(mfrow=c(3,6), mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0), tck=-0.02, cex.axis=1.2)
for (lats in 1:dim(rad2.diff.mon)[2]){ 
  plot(rad2.diff.mon[,lats,1,1], type='n',ylim=range(rad2.diff.mon[,,1,1]))
  grid(col='gray80')
  legend('topright',legend=paste('lattitude =',lat0[lats],sep=''),inset=0,bty='n') #,bg='gray80'
  points(rad2.diff.mon[,lats,1,1],pch=16, cex=.7)
}
dev.off()

SS <- readMat('../../codes/matLM/SS.mat')
labCreal <- as.numeric(SS$labsC)

# get figure rdm
postscript('realdatamap.eps', pointsize=30,width=16,height=12,horizontal=F,paper='special')
Yreal <- matrix(NA,T,N)
varid <- 1
ran <- range(rad2.diff.mon[inds, , varid, ])
dat <- matrix(NA, tnf, nx*ny)
scal1 <- 1.2; scal2 <- .00025; scal3 <- max(rad2.diff.mon[inds, , 1, ])
k <- 1
palette(rainbow(max(labCreal)))
par(mar=c(2.4,2.4,0,0)+.05,mgp = c(1.3, .2, 0),tck=-0.02,cex.axis=1.2,cex.lab=1.2)
plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude',main=''); #grid(col='gray70')
for(i in 1:nx){
  for(j in 1:ny){
    Yreal[,k] <- y <- rad2.diff.mon[inds, j, varid, i]
    y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
    x1 <- c(1:length(y1))*scal2+month[i]
    lines(y1 ~ x1, lwd=1, col=max(labCreal)+1-labCreal[k]) #col='gray20',
    #bp <- round(Box.test(y, type="Ljung-Box", lag=1)$p.value,3)
    #text(x1[1],y1[1]-.01,labels=as.character(bp),col='gray20',cex=.8)
    k <- k+1
  }
}
dev.off()

# get figure rdm2
labCreal2 <- labmat[,1]
table(labCreal2,labCreal)
postscript('realdatamap2.eps', pointsize=30,width=16,height=12,horizontal=F,paper='special')
Yreal <- matrix(NA,T,N)
varid <- 1
ran <- range(rad2.diff.mon[inds, , varid, ])
dat <- matrix(NA, tnf, nx*ny)
scal1 <- 1.2; scal2 <- .00025; scal3 <- max(rad2.diff.mon[inds, , 1, ])
k <- 1
palette(rainbow(max(labCreal2)))
par(mar=c(2.4,2.4,0,0)+.05,mgp = c(1.3, .2, 0),tck=-0.02,cex.axis=1.2,cex.lab=1.2)
plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude',main=''); #grid(col='gray70')
for(i in 1:nx){
  for(j in 1:ny){
    Yreal[,k] <- y <- rad2.diff.mon[inds, j, varid, i]
    y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
    x1 <- c(1:length(y1))*scal2+month[i]
    lines(y1 ~ x1, lwd=1, col=max(labCreal2)+1-labCreal2[k]) #col='gray20',
    #bp <- round(Box.test(y, type="Ljung-Box", lag=1)$p.value,3)
    #text(x1[1],y1[1]-.01,labels=as.character(bp),col='gray20',cex=.8)
    k <- k+1
  }
}
dev.off()

# get figure rdmr
Yreal <- matrix(NA,T,N)
postscript('realdatamapraw.eps', pointsize=30,width=16,height=12,horizontal=F,paper='special')
varid <- 1
ran <- range(rad2.diff.mon[inds, , varid, ])
dat <- matrix(NA, tnf, nx*ny)
scal1 <- 1.2; scal2 <- .00025; scal3 <- max(rad2.diff.mon[inds, , 1, ])
k <- 1
palette(rainbow(d))
par(mar=c(2.4,2.4,0,0)+.05,mgp = c(1.3, .2, 0),tck=-0.02,cex.axis=1.2,cex.lab=1.2)
plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude',main=''); #grid(col='gray70')
for(i in 1:nx){
  for(j in 1:ny){
    Yreal[,k] <- y <- rad2.diff.mon[inds, j, varid, i]
    y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
    x1 <- c(1:length(y1))*scal2+month[i]
    lines(y1 ~ x1, lwd=1) #col='gray20',
    #bp <- round(Box.test(y, type="Ljung-Box", lag=1)$p.value,3)
    #text(x1[1],y1[1]-.01,labels=as.character(bp),col='gray20',cex=.8)
    k <- k+1
  }
}
dev.off()

# get figure rdc: with mean of sites in cluster
postscript(file='realdatabycluster.eps',pointsize=25,width=16,height=8,horizontal=F,paper='special')
palette(rainbow(max(labCreal)))
par(mfrow=c(2,2), mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0), tck=-0.02, cex.axis=1.2)
for(r in 1:max(labCreal)){
  indr <- which(labCreal==r)
  matplot(Yreal[,indr], type='n', ylim=range(Yreal)); for(s in indr)  lines(1:T, Yreal[,s], col='gray70')
  # betar <- as.numeric(simus$truePara[[6]])
  # lines(1:T, waveT%*%betar[(r-1)*T+(1:T)], lwd=3, col='red')
  lines(1:T, rowMeans(Yreal[,indr]), lwd=3, col=max(labCreal)+1-r)
}
dev.off()

# get figure rdc2: with mean and confidence band
bands <- readMat('../../codes/curvebands.mat')$mat
postscript(file='realdatabycluster2.eps',pointsize=25,width=16,height=8,horizontal=F,paper='special')
palette(rainbow(max(labCreal)))
par(mfrow=c(2,2), mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0), tck=-0.02, cex.axis=1.2)
for(r in 1:max(labCreal)){
  indr <- which(labCreal==r)
  matplot(Yreal[,indr], type='n', ylim=range(Yreal)); for(s in indr)  lines(1:T, Yreal[,s], col='gray70')
  lines(1:T, bands[(r-1)*T+(1:T),1], lwd=3, col=max(labCreal)+1-r)
  lines(1:T, bands[(r-1)*T+(1:T),2], lwd=2, col='black', lty=1)
  lines(1:T, bands[(r-1)*T+(1:T),3], lwd=2, col='black', lty=1)
}
dev.off()


# get Figure rdss
postscript(file='rdss.eps',pointsize=25,width=10,height=8,horizontal=F,paper='special')
par(mar=c(2.4,2.4,0,0)+.05,mgp = c(1.3, .2, 0), tck=-0.01, cex.axis=1.2, cex.lab=1.2)
plot(SS$Sensitivity, SS$Specificity,type='n',xlim=c(0,1),ylim=c(0,1),xlab='Sensitivity', ylab='Specificity')  
grid(col='gray50'); points(SS$Sensitivity, SS$Specificity, pch=1, cex=.9, xlab='Sensitivity', ylab='Specificity', col='gray30')
abline(h=1, lwd=3, lty=2, col='red');  abline(v=1, lwd=3, lty=2, col='red')
dev.off()


# ---------------------------------------------------- SIMULATION
# Figure sdm
postscript('simulatedatamap.eps', pointsize=30,width=16,height=12,horizontal=F,paper='special')
ran <- range(Y0)
scal1 <- 1.5; scal2 <- .00026; scal3 <- max(Y0)
palette(rainbow(d))
wds <- rep(1,N); wds[cts] <- 3
par(mar=c(2.4,2.4,0,0)+.05,mgp = c(1.3, .2, 0),tck=-0.02,cex.axis=1.2,cex.lab=1.2)
plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude'); #grid(col='gray70')
s <- 1
for(i in 1:nx){
  for(j in 1:ny){
    y <- Y0[,s]
    y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
    x1 <- c(1:length(y1))*scal2+month[i]
    lines(y1 ~ x1, lwd=wds[s], col=labs[s]) #col='gray20',
    s <- s+1
  }
}
dev.off()

SS7 <- readMat('../../codes/matLM/simu/SS7.mat')
labCsimu <- as.integer(SS7$labsC)

# Figure sdm2
(tab <- table(labCsimu,labs))
vec <- integer(d); for(r in 1:d) vec[r] <- which(tab[r,]==max(tab[r,]))

postscript('simulatedatamapMODEL.eps', pointsize=30,width=16,height=12,horizontal=F,paper='special')
ran <- range(Y0)
scal1 <- 1.5; scal2 <- .00026; scal3 <- max(Y0)
palette(rainbow(d))
par(mar=c(2.4,2.4,0,0)+.05,mgp = c(1.3, .2, 0),tck=-0.02,cex.axis=1.2,cex.lab=1.2)
plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude'); #grid(col='gray70')
s <- 1
for(i in 1:nx){
  for(j in 1:ny){
    y <- Y0[,s]
    y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
    x1 <- c(1:length(y1))*scal2+month[i]
    lines(y1 ~ x1, lwd=1, col=vec[labCsimu[s]]) #col='gray20',
    s <- s+1
  }
}
dev.off()


# get Figure sdss
postscript(file='sdss.eps',pointsize=25,width=10,height=8,horizontal=F,paper='special')
par(mar=c(2.4,2.4,0,0)+.05,mgp = c(1.3, .2, 0), tck=-0.01, cex.axis=1.2, cex.lab=1.2)
plot(SS7$SensitivityC, SS7$SpecificityC,type='n',xlim=c(0,1),ylim=c(0,1),xlab='Sensitivity', ylab='Specificity')  
grid(col='gray50'); points(SS7$SensitivityC, SS7$SpecificityC, pch=1, cex=.9, xlab='Sensitivity', ylab='Specificity', col='gray30')
abline(h=1, lwd=3, lty=2, col='red');  abline(v=1, lwd=3, lty=2, col='red')
dev.off()


# get Figure sdss2
postscript(file='sdss2.eps',pointsize=25,width=10,height=8,horizontal=F,paper='special')
par(mar=c(2.4,2.4,0,0)+.05,mgp = c(1.3, .2, 0), tck=-0.01, cex.axis=1.2, cex.lab=1.2)
plot(SS7$Sensitivity, SS7$Specificity,type='n',xlim=c(0,1),ylim=c(0,1),xlab='Sensitivity', ylab='Specificity')  
grid(col='gray50'); points(SS7$Sensitivity, SS7$Specificity, pch=1, cex=.9, xlab='Sensitivity', ylab='Specificity', col='gray30')
abline(h=1, lwd=3, lty=2, col='red');  abline(v=1, lwd=3, lty=2, col='red')
points(SS7$SensitivityCT, SS7$SpecificityCT, pch=16, col='red',cex=1.6)
dev.off()

# get Figure sdc
postscript(file='simulatedatabycluster.eps',pointsize=25,width=16,height=8,horizontal=F,paper='special')
par(mfrow=c(2,4), mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0), tck=-0.02, cex.axis=1.2)
for(r in 1:d){
  indr <- which(labs==r)
  matplot(Y0[,indr], type='n', ylim=range(Y0)); for(s in indr)  lines(1:T, Y0[,s], col='gray70')
  betar <- as.numeric(simus$truePara[[6]])
  lines(1:T, waveT%*%betar[(r-1)*T+(1:T)], lwd=3, col='red')
}
dev.off()


# get figure rdc: with mean of sites in cluster
postscript(file='simudatabycluster.eps',pointsize=25,width=16,height=8,horizontal=F,paper='special')
palette(rainbow(max(labCsimu)))
par(mfrow=c(2,3), mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0), tck=-0.02, cex.axis=1.2)
for(r in 1:max(labCsimu)){
  indr <- which(labCsimu==r)
  matplot(Y0[,indr], type='n', ylim=range(Yreal)); for(s in indr)  lines(1:T, Y0[,s], col='gray70')
  lines(1:T, rowMeans(Y0[,indr]), lwd=3, col=vec[r])
}
dev.off()

# get figure rdc2: with mean and confidence band
bands7 <- readMat('../../codes/matLM/simu/curvebands7.mat')$mat
postscript(file='simudatabycluster2.eps',pointsize=25,width=16,height=8,horizontal=F,paper='special')
palette(rainbow(max(labCsimu)))
par(mfrow=c(2,3), mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0), tck=-0.02, cex.axis=1.2)
for(r in 1:max(labCsimu)){
  indr <- which(labCsimu==r)
  matplot(Y0[,indr], type='n', ylim=range(Y0)); for(s in indr)  lines(1:T, Y0[,s], col='gray70')
  lines(1:T, bands7[(r-1)*T+(1:T),1], lwd=3, col=vec[r])
  lines(1:T, bands7[(r-1)*T+(1:T),2], lwd=2, col='black', lty=1)
  lines(1:T, bands7[(r-1)*T+(1:T),3], lwd=2, col='black', lty=1)
}
dev.off()


# --------------------------------------------------------------------------------------
# get other figures
# G plot

Gs <- readMat('../../codes/G.mat')

#pdf(file='Growth_x1_cl.pdf',pointsize=25,width=8,height=6)
postscript(file='Gx1_cl.eps',pointsize=25,width=8,height=4,horizontal=F,paper='special')
labsx1 <- Gs$x1[[3]]
palette(rainbow(max(labsx1)))
par(mfrow=c(2,3), mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0), tck=-0.02, cex.axis=1.2)
for (r in 1:as.numeric(Gs$x1[[1]])) { 
  yobs <- Y[,labsx1==r]
  yhat <- Gs$X%*%Gs$x1[[7]][(r-1)*T+(1:T)]
  matplot(yobs,type='n', ylim=range(cbind(yobs, yhat,Gs$X%*%Gs$x[[7]][(r-1)*T+(1:T)]),na.rm=T)); for(i in 1:ncol(yobs)) points(1:T,yobs[,i],pch=16,col='gray50'); lines(1:T,yhat,lwd=2, col=r) 
}   
dev.off()

#pdf(file='Growth_x1_map.pdf',pointsize=25,width=8,height=6)
postscript('Gx1_map.eps', pointsize=30,width=8,height=6,horizontal=F,paper='special')
k <- 1
wds <- rep(1,N); wds[Gs$x1[[2]]] <- 3
palette(rainbow(max(labsx1)))
par(mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0),cex.axis=1);
plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude',main=''); #grid(col='gray70')
for(i in 1:nx){
  for(j in 1:ny){
    y <- rad2.diff.mon[inds, j, varid, i]
    y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
    x1 <- c(1:length(y1))*scal2+month[i]
    lines(y1 ~ x1, lwd=wds[k], col=labsx1[k])
    k <- k+1
  }
}
dev.off()

#pdf(file='Growth_x_cl.pdf',pointsize=25,width=8,height=6)
postscript(file='Gx_cl.eps',pointsize=25,width=8,height=4,horizontal=F,paper='special')
labsx <- Gs$x[[3]]
palette(rainbow(max(labsx)))
par(mfrow=c(2,3), mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0), tck=-0.02, cex.axis=1.2)
for (r in 1:as.numeric(Gs$x[[1]])) { 
  yobs <- Y[,labsx==r]
  yhat <- Gs$X%*%Gs$x[[7]][(r-1)*T+(1:T)] 
  matplot(yobs,type='n',ylim=range(cbind(yobs, yhat,Gs$X%*%Gs$x1[[7]][(r-1)*T+(1:T)]),na.rm=T)); for(i in 1:ncol(yobs)) points(1:T,yobs[,i],pch=16,col='gray50'); lines(1:T,yhat,lwd=2, col=r) 
}   
dev.off()

#pdf(file='Growth_x_map.pdf',pointsize=25,width=8,height=6)
postscript('Gx_map.eps', pointsize=30,width=8,height=6,horizontal=F,paper='special')
wds <- rep(1,N); wds[Gs$x[[2]]] <- 3
k <- 1
palette(rainbow(max(labsx)))
par(mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0),cex.axis=1);
plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude',main=''); #grid(col='gray70')
for(i in 1:nx){
  for(j in 1:ny){
    y <- rad2.diff.mon[inds, j, varid, i]
    y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
    x1 <- c(1:length(y1))*scal2+month[i]
    lines(y1 ~ x1, col=labsx[k],lwd=wds[k])
    k <- k+1
  }
}
#dev.off()
dev.off()






Gs <- readMat('../../codes/M.mat')

#pdf(file='Growth_x1_cl.pdf',pointsize=25,width=8,height=6)
postscript(file='Mx1_cl.eps',pointsize=25,width=8,height=4,horizontal=F,paper='special')
labsx1 <- Gs$x1[[3]]
palette(rainbow(max(labsx1)))
par(mfrow=c(2,3), mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0), tck=-0.02, cex.axis=1.2)
for (r in 1:as.numeric(Gs$x1[[1]])) { 
  yobs <- Y[,labsx1==r]
  yhat <- Gs$X%*%Gs$x1[[7]][(r-1)*T+(1:T)]
  matplot(yobs,type='n', ylim=range(cbind(yobs, yhat,Gs$X%*%Gs$x[[7]][(r-1)*T+(1:T)]),na.rm=T)); for(i in 1:ncol(yobs)) points(1:T,yobs[,i],pch=16,col='gray50'); lines(1:T,yhat,lwd=2, col=r) 
}   
dev.off()

#pdf(file='Growth_x1_map.pdf',pointsize=25,width=8,height=6)
postscript('Mx1_map.eps', pointsize=30,width=8,height=6,horizontal=F,paper='special')
k <- 1
wds <- rep(1,N); wds[Gs$x1[[2]]] <- 3
palette(rainbow(max(labsx1)))
par(mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0),cex.axis=1);
plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude',main=''); #grid(col='gray70')
for(i in 1:nx){
  for(j in 1:ny){
    y <- rad2.diff.mon[inds, j, varid, i]
    y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
    x1 <- c(1:length(y1))*scal2+month[i]
    lines(y1 ~ x1, lwd=wds[k], col=labsx1[k])
    k <- k+1
  }
}
dev.off()

#pdf(file='Growth_x_cl.pdf',pointsize=25,width=8,height=6)
postscript(file='Mx_cl.eps',pointsize=25,width=8,height=4,horizontal=F,paper='special')
labsx <- Gs$x[[3]]
palette(rainbow(max(labsx)))
par(mfrow=c(2,3), mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0), tck=-0.02, cex.axis=1.2)
for (r in 1:as.numeric(Gs$x[[1]])) { 
  yobs <- Y[,labsx==r]
  yhat <- Gs$X%*%Gs$x[[7]][(r-1)*T+(1:T)] 
  matplot(yobs,type='n',ylim=range(cbind(yobs, yhat,Gs$X%*%Gs$x1[[7]][(r-1)*T+(1:T)]),na.rm=T)); for(i in 1:ncol(yobs)) points(1:T,yobs[,i],pch=16,col='gray50'); lines(1:T,yhat,lwd=2, col=r) 
}   
dev.off()

#pdf(file='Growth_x_map.pdf',pointsize=25,width=8,height=6)
postscript('Mx_map.eps', pointsize=30,width=8,height=6,horizontal=F,paper='special')
wds <- rep(1,N); wds[Gs$x[[2]]] <- 3
k <- 1
palette(rainbow(max(labsx)))
par(mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0),cex.axis=1);
plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude',main=''); #grid(col='gray70')
for(i in 1:nx){
  for(j in 1:ny){
    y <- rad2.diff.mon[inds, j, varid, i]
    y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
    x1 <- c(1:length(y1))*scal2+month[i]
    lines(y1 ~ x1, col=labsx[k],lwd=wds[k])
    k <- k+1
  }
}
#dev.off()
dev.off()






Gs <- readMat('../../codes/S.mat')

#pdf(file='Growth_x1_cl.pdf',pointsize=25,width=8,height=6)
postscript(file='Sx1_cl.eps',pointsize=25,width=8,height=4,horizontal=F,paper='special')
labsx1 <- Gs$x1[[3]]
palette(rainbow(max(labsx1)))
par(mfrow=c(2,3), mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0), tck=-0.02, cex.axis=1.2)
for (r in 1:as.numeric(Gs$x1[[1]])) { 
  yobs <- Y[,labsx1==r]
  yhat <- Gs$X%*%Gs$x1[[7]][(r-1)*T+(1:T)]
  matplot(yobs,type='n', ylim=range(cbind(yobs, yhat,Gs$X%*%Gs$x[[7]][(r-1)*T+(1:T)]),na.rm=T)); for(i in 1:ncol(yobs)) points(1:T,yobs[,i],pch=16,col='gray50'); lines(1:T,yhat,lwd=2, col=r) 
}   
dev.off()

#pdf(file='Growth_x1_map.pdf',pointsize=25,width=8,height=6)
postscript('Sx1_map.eps', pointsize=30,width=8,height=6,horizontal=F,paper='special')
k <- 1
wds <- rep(1,N); wds[Gs$x1[[2]]] <- 3
palette(rainbow(max(labsx1)))
par(mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0),cex.axis=1);
plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude',main=''); #grid(col='gray70')
for(i in 1:nx){
  for(j in 1:ny){
    y <- rad2.diff.mon[inds, j, varid, i]
    y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
    x1 <- c(1:length(y1))*scal2+month[i]
    lines(y1 ~ x1, lwd=wds[k], col=labsx1[k])
    k <- k+1
  }
}
dev.off()

#pdf(file='Growth_x_cl.pdf',pointsize=25,width=8,height=6)
postscript(file='Sx_cl.eps',pointsize=25,width=8,height=4,horizontal=F,paper='special')
labsx <- Gs$x[[3]]
palette(rainbow(max(labsx)))
par(mfrow=c(2,3), mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0), tck=-0.02, cex.axis=1.2)
for (r in 1:as.numeric(Gs$x[[1]])) { 
  yobs <- Y[,labsx==r]
  yhat <- Gs$X%*%Gs$x[[7]][(r-1)*T+(1:T)] 
  matplot(yobs,type='n',ylim=range(cbind(yobs, yhat,Gs$X%*%Gs$x1[[7]][(r-1)*T+(1:T)]),na.rm=T)); for(i in 1:ncol(yobs)) points(1:T,yobs[,i],pch=16,col='gray50'); lines(1:T,yhat,lwd=2, col=r) 
}   
dev.off()

#pdf(file='Growth_x_map.pdf',pointsize=25,width=8,height=6)
postscript('Sx_map.eps', pointsize=30,width=8,height=6,horizontal=F,paper='special')
wds <- rep(1,N); wds[Gs$x[[2]]] <- 3
k <- 1
palette(rainbow(max(labsx)))
par(mar=c(1.4,1.4,0,0)+.05,mgp = c(3, .2, 0),cex.axis=1);
plot(1, type='n', ylim=range(lat), xlim=range(month),xlab='month',ylab='latitude',main=''); #grid(col='gray70')
for(i in 1:nx){
  for(j in 1:ny){
    y <- rad2.diff.mon[inds, j, varid, i]
    y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j]
    x1 <- c(1:length(y1))*scal2+month[i]
    lines(y1 ~ x1, col=labsx[k],lwd=wds[k])
    k <- k+1
  }
}
#dev.off()
dev.off()



