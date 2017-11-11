
# 2017-10-26: useful functions

#++++++++++++++++++++++++++ setting for global variables
nx <- length(month <- c(1:24))
ny <- length(lat <- seq(-85,85, by=10))

lat2lab <- function(lat0, padding=FALSE){
  # value of latitude to label with degree and S/N indicator
  inds <- which(grepl('-', lat0))
  sn <- rep('N',length(lat0)); sn[inds] <- 'S'
  if(padding) sn <- paste(sn, "              ")
  tmp <- gsub('-','',lat0)
  lab <- expression(); for(i in 1:length(lat0)) lab[i] <- as.expression(substitute(paste(a1,degree,a2),list(a1=tmp[i], a2=sn[i])))
  return(lab)
}


sortLab <- function(lab){
  # sort cluster labels so that 1 = South Pole, max = North Pole
  lab <- factor(lab,levels=1:max(lab))
  lats <- rep(lat, nx); minLat <- aggregate(lats, list(lab), 'min')
  x <- minLat$x
  w <- aggregate(lats, list(lab), function(x) sum(x==min(x)))$x / 1e8
  x <- x - w # the more weight (frequency of the minimal values), the smaller x to get rid of ties
  lab1 <- as.integer(factor(lab, levels=order(x)))
  return(lab1)
}


plotMapLab <- function(labs, Y, varids=1, saveFig=FALSE, fnam=NULL){
  # labs is N-vector for cluster labels
  # Y is T by N data matrix
  
  if(saveFig){
    # postscript(fnam, width=6.8, height=5.4,  pointsize=11, horizontal=F,paper='special')
    #jpeg(fnam,quality=90, height=5.4,width=6.8, units='in', pointsize=11, res=300)
    #pdf(fnam, width = 6.2, height = 5.4, pointsize=11,paper='special')
    png(fnam, height=5,width=6.8, units='in', pointsize=12, res=600)
  }
  
  Ynams <- c('Spectral radiance (control run)', 'Near-surface and skin temperature', 'Surface-200-hPa temperature', '200-10-hPa temperature', 'Surface-500-hPa water vapor', '500-200-hPa water vapor', 'Low-level cloud-top height', 'Midlevel cloud-top height', 'High-level cloud-top height', 'Low-level cloud fraction', 'Midlevel cloud fraction', 'High-level cloud fraction', 'Thin ice cloud optical thickness', 'Ice cloud particle size', 'Water cloud optical thickness', 'Water cloud particle size')
  
  # Ynams <-  c('Control run', 'Near-surface&skin temperature', 'Surface-200-hPa temperature', '200-10-hPa temperature', 'Surface-500-hPa water vapor', '500-200-hPa water vapor', 'Low-level cloud-top height', 'Midlevel cloud-top height', 'High-level cloud-top height', 'Low-level cloud fraction', 'Midlevel cloud fraction', 'High-level cloud fraction', 'Thin ice cloud optical thickness', 'Ice cloud particle size', 'Water cloud optical thickness', 'Water cloud particle size')
  
  # require(RColorBrewer)
  # cols <- c(brewer.pal(9,'Set1'),brewer.pal(12,'Paired'),brewer.pal(8,'Set2'))
  # cols[6] <- '#fee08b'
  # cols <- cols[c(4,6,5,1,3,2)]
  # cols <- c("#984EA3", "#fee08b", "#FF7F00", "#E41A1C", "#4DAF4A", "#377EB8")
  # cols <- c("#984EA3", "#a65628", "#FF7F00", "#E41A1C", "#4DAF4A", "#377EB8")
  cols <- c('#984ea3','#a65628','#ff7f00','#e41a1c','#4daf4a','#377eb8','#4d4d4d','#c51b7d','#016c59','#a50f15','#feb24c','#253494')
  # plot(1:length(cols), col=cols,pch=16)
  
  for (varid in varids) {
    if(varid==1) ran <- range(Y) else ran <- range(Y[varid-1])
    
    scal1 <- 11; scal2 <- .007; 
    
    
    ### when the X were multiplied by 1e4 --- prior to 2017-10-30 data processing
    # scal1 <- 80; scal2 <- .006; 
    # if(varid%in%c(10)) scal1 <- 15e3; scal2 <- .006; 
    # if(varid%in%c(2,3,12)) scal1 <- 35e3; scal2 <- .006; 
    # if(varid%in%c(9)) scal1 <- 42e3; scal2 <- .006; 
    # if(varid%in%c(4,8)) scal1 <- 75e3; scal2 <- .006; 
    # if(varid%in%c(13)) scal1 <- 3e5; scal2 <- .006;
    # if(varid%in%c(14)) scal1 <- 6e5; scal2 <- .006;
    # if(varid%in%c(6)) scal1 <- 8e5; scal2 <- .006; 
    # if(varid%in%c(11,15,16)) scal1 <- 1e6; scal2 <- .006; 
    # if(varid%in%c(7)) scal1 <- 5e6; scal2 <- .006; 
    # if(varid%in%c(5)) scal1 <- 4e7; scal2 <- .006; 
    
    scal3 <- max(Y)
    # par(mar=c(1.9,2,1.2,.2)+0,mgp=c(1,.2,0), tck=-0.01, cex.axis=.85, cex.lab=1, cex.main=.8, xaxt='n', yaxt='n')
    m0 <- matrix(1:2, nrow=2, ncol=1, byrow=TRUE)
    layout(mat=m0, heights=c(.94,.06))
    par(oma=c(0,0,0,0), mar=c(1.9,2,1.2,.2)+0,mgp=c(1,.2,0), tck=-0.01, cex.axis=.85, cex.lab=1, cex.main=.8, xaxt='n', yaxt='n')
    
    plot(1, type='n', 
         ylim=range(lat)+c(1,.5), 
         xlim=range(month)+c(.5,0),
         xlab='',#'month',
         ylab='', #'latitude',
         main=Ynams[varid] 
    )
    abline(v=month,col='gray70', lty=3)
    abline(h=lat,col='gray70', lty=3)
    k <- 1
    for(i in 1:nx){
      for(j in 1:ny){
        y <- Y[,k]
        y1 <- (y - y[1])*scal1*scal3/(ran[2]-ran[1])+lat[j] - 2.5
        x1 <- c(1:length(y1))*scal2+month[i] - .3
        lines(y1 ~ x1, lwd=1, col=cols[labs[k]]) #,d+1-labmat[k,varid]
        k <- k+1
      }
    }
    par(xaxt='s', yaxt='s')
    par(las=1,mgp=c(1.2,.2,0), tck=-0.005); 
    mytime <- paste(rep(c(2003,2004), each=12), rep(1:12,2),sep='/')
    # axis(side=1,at=1:24,labels=mytime) #seq(1,24,by=1)
    axis(side=1,at=1:24,labels=NA) 
    text(c(1:24)+.5,rep(-90,24), labels = mytime ,srt = 40,adj = c(1.1,1.1),xpd = TRUE, cex=.7)
    par(las=1,mgp=c(1.2,.2,0), tck=-0.005); 
    axis(side=2,at=lat,labels=lat2lab(lat))
    
    # check
    a <- as.numeric(readMat('checkFindSP.mat')$labs)
    text(x=rep(1:24, each=18), y=rep(seq(-85,85, by=10),24), labels=a)
    
    labNam <- 1:max(labs)
    par(mar=c(0,0,0,0))
    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    legend('center',legend=labNam, bty='n', fill=NA, border=NA, text.col='black',inset=0, ncol=length(labNam), cex=.8, pt.cex=15, pch='.',col=cols)  #, title='Cluster'
  }
  
  if(saveFig) dev.off()
  
}



# not run


