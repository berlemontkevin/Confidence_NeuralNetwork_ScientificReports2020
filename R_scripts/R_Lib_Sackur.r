library('MASS')
library('pracma') # needed for the geometry in ROC curves.
#library('psyphy')
#library('sfsmisc')
#library('Hmisc')
## library('gplots')
## library('gdata')
## library('caTools')
## library('gtools')
# library(rjags)
# load.module("wiener")
#source('~/science/0librairies/chrplr.r')

# Pour réinstaller la suite des paquets que tu utilises. Lancer une session R avec les droits d'administrateur, puis
# install.packages(c("fields","ggplot2", "gtools",  "proto", "psy", "psych", "psychometric", "sfsmisc", "spam", "bootstrap", "boot", "caTools", "gsubfn", "subplex", "smacof"))


dmatplot <- function(x,y,errorbars,dep,errlen=.05,ylim,...) {

  n=ncol(y)
  b=rep(1,length(x))

  if (missing(dep)) { dep=min(diff(sort(x)))/10 }
  s=seq(0,by=dep,length=n)

  x2=x + t(s) %x% b

  if (missing (errorbars)) errorbars=0
  if (missing (ylim)) ylim= c (min (y-errorbars), max (y+errorbars))
  
  matplot(x2, y, ylim=ylim, col=1, axes=F, frame.plot=T, type="b",...)

  if (length(errorbars)>1)
  {
    x0=as.vector(x2)
    y0=as.vector(y)
    err=as.vector(errorbars)
    arrows(x0,y0-err,x0,y0+err,code=3,angle=90,length=errlen)
  }
  structure (list(ylim=ylim))
}


rm.main.effect <- function (x,effect,FUN=mean) {
  m=tapply(x,effect,FUN)
  x-m[as.character(effect)]
}
numerize <- function (x){
  return (as.numeric (as.character (x)))}

se <- function(x) {
    nx <- length(x)
    if (nx < 2) 
        stop("not enough x observations")
    sqrt(var(x, na.rm=T)/nx)
}

sat.dens <- function (rt, acc, norm=T, out=.05, kernel='rectangular') {
  q <- quantile (rt, probs=c (out, 1-out))
  rtn <- rt[rt>q[1]&rt<q[2]]
  acc <- acc[rt>q[1]&rt<q[2]]
  rt <- rtn
  de <- density(rt[acc==0], kernel=kernel)
  dc <- density(rt[acc==1], kernel=kernel)
  ref <- findInterval(de$x, dc$x)
  lmax <- min ((1:512)[ref==max (ref)])
  ref <- ref[1:lmax]
  ref <- ref[ref!=0]
  dcx <- dc$x[ref]
  dex <- de$x[1:length(ref)]
  if (norm==T) {
    dc$y <- dc$y/max (dc$y)
    de$y <- de$y/max (de$y)
  }
  tacc <- 1-de$y[1:length(ref)]/(dc$y[ref]+de$y[1:length(ref)])
  structure (list (x=dcx, tacc=tacc))
}


sat.qt <- function (rt, acc, n) {
  ref <- gl(n, ceiling(length(rt)/n), length (rt))
  ref <- as.factor(ref[rank(abs(rt))])
  tacc <- tapply(acc, list (ref=ref), mean)
  x <- as.vector(quantile (rt, probs=seq (0, 1, by=(1/n))))
  x <- (diff (x)/2)+ x[1:(length(x)-1)]
  structure (list (x=x, tacc=tacc))
}


sat.abs <- function (rt, acc, vec) {
#   Cette fonction calcule les taux d'erreur _A L'intérieur_ des
#   intervales donnés par vec. Et exclu le reste, comme les deux
#   premières lignes le montrent.
  mi <- vec[1];ma <- vec[length(vec)]
  rtt <- abs (rt[rt<=ma&rt>mi])
  acc <- acc[rt<=ma&rt>mi]
  ref <- findInterval (vec, sort (rtt))
  levels <- rep (1:(length(ref)-1), diff(ref))
  u <- as.factor(levels[rank(rtt)])
  tacc <- tapply(acc, list (u=u), mean)
  x <- (diff (vec)/2)+ vec[1:(length(vec)-1)]
  structure (list (x=x, tacc=tacc))
}  


### Version de Christophe Pallier

sat <- function (rt, acc,n)
  {
     breaks = quantile(rt,probs=seq(0,1,1/n))
     centers = (breaks[-1] + breaks[-length(breaks)])/2
     accs = tapply(acc, cut(rt, breaks ), mean)
     list(centers=centers, accuracies=accs)
  }

bin.count <- function (var, resp, n, responses=c(0,1)) {
  ref <- gl(n, ceiling(length(var)/n), length (var))
  ref <- as.factor(ref[rank(abs(var))])
  t0 <- tapply(resp==responses[1], list (ref=ref), sum)
  t1 <- tapply(resp==responses[2], list (ref=ref), sum)
  x <- as.vector(quantile (var, probs=seq (0, 1, by=(1/n))))
  x <- (diff (x)/2)+ x[1:(length(x)-1)]
  t <- data.frame(x=x,t0=t0,t1=t1)
  t$phat <- t$t1/(t$t1+t$t0)
  structure (t)
}

bin.count.fixed <- function (var, resp, n, responses=c(0,1)) {
  breaks <- seq(min(var), max(var), by=(max(var)-min(var))/n)
  centers <- (breaks[-1] + breaks[-length(breaks)])/2
  t0 <- tapply(resp==responses[1], cut(var, breaks), sum)
  t1 <- tapply(resp==responses[2], cut(var, breaks), sum)
  t <- data.frame(x=centers,t0=t0,t1=t1)
  structure (t)
}



nbmodes <- function (x, h)
# returns how many modes there are in the kernel density estimate of
# vector x, with window width h.
  {
    modes <- density (x, bw=h)
    modes <- diff (diff (modes$y) / abs (diff (modes$y)))
    return (sum(modes==-2))
  }

hcrit <- function (x, n, e=.0001)
# Returns the minimal window width such that the kernel density estimate of x
# has n modes. 
{
  minw <- min (abs (diff (x)))
  maxw <- (max (x) - min (x))/2
  winN <- maxw
  winN1 <- minw
  while (abs(winN-winN1)>e)
    {
      modes <- nbmodes (x, winN)
      if (is.na(modes)) print(paste(winN, winN1, maxw, minw))
      winN1 <- winN
      if (modes > n)
        {
          minw <- winN
          winN <- (winN + maxw)/2
        }
      else
        {
          maxw <- winN
          winN <- (winN + minw)/2
        }
    }
return (winN)
}


silvers <- function (x, m, nboot=200)
# Silverman's significance for the null that x has at
# most m modes.
{
  h0 <- hcrit (x, m)
  n <- 0
  l <- length(x)
  c <- sqrt ( 1+ (h0^2/var (x)))
  for (i in 1:nboot) {
    boot <- (sample(x, replace=T) + rnorm (l, mean=0, sd=h0))/c
    nb <- nbmodes (boot, h0)
    if (nb > m) {
      n <- n+1
    }
  }
  return (n/nboot)
}

pdip <- function (dipstat, size){
#  data (qDiptab)
  tSizes <- as.numeric (names (qDiptab [,1]))
  rc <- diff (rank (c (tSizes, size)))
  rankSize <- ((1:(length (tSizes)-1))[rc==2])[1]
  prop <- (size-tSizes[rankSize])/(tSizes[rankSize+1]-tSizes[rankSize])
  z <- qDiptab[rankSize,]+(prop*(qDiptab[rankSize+1,]-qDiptab[rankSize,]))
  t <- as.numeric (names (qDiptab [1,]))
  fonction <- approxfun(z, t)
  pval <- fonction(dipstat)
  return(pval)
}

apdip <- function (x){
  pdip (dip (x), length (x))
}


# Chondrite meteor data
chond <- c(20.77, 22.56, 22.71, 22.99, 26.39, 27.08, 27.32, 27.33, 27.57,
       27.81, 28.69, 29.36, 30.25, 31.89, 32.88, 33.23, 33.28, 33.40, 33.52,
       33.83, 33.95, 34.82)
 
locmodes <- function (x, n)

  {
    h <- hcrit(x, n)
    d <- density (x, bw=h+.0001)
    modes <- diff (diff (d$y) / abs (diff (d$y)))
    x <- d$x[modes==-2]
    y <- d$y[modes==-2]
    return (list(pos=x, freq=y))
  }

siegel <- function (x,y) {
  library(ctest)
  
  f <- function (n,l,cut) {
    x <- n-cut
    y <- cut-n
    m <- ((abs(x)/x+1)/2)*(2*(l-n)+2-((l-n)%%2)) +
      ((abs(y)/y+1)/2)*(2*n-(n%%2))
    m
  }

  z <- c (x, y)
  znames <- c (rep("x", length (x)), rep("y", length (y)))
  zranks <- rank (z)
  l <- length (z)
  cut <- (l-(l%%2))/2
  zreord <- f(zranks, l, cut)
  zreord[is.nan(zreord)] <- 2*cut-(cut%%2)
  d <- data.frame(znames, zranks, zreord)
  wilcox.test(d$zreord[d$znames=="x"],d$zreord[d$znames=="y"])
}

moses <- function (x,y,s) {
  x <- sample (x)
  y <- sample (y)
  xl <- length (x)%/%s
  yl <- length (y)%/%s
  x <- x[1:(xl*s)]
  y <- y[1:(yl*s)]
  xf <- gl (xl, s)
  yf <- gl (yl, s)
  xr <- aggregate (x, list (xf=xf), var)
  yr <- aggregate (y, list (yf=yf), var)
  wilcox.test (xr$x, yr$x)
}


#===============================================================================#
# Calcule un vecteur qui vaut T lorsqu'il y a répétition à intervalle n         #
# dans le vecteur x. Sauf si les valeurs dans le vecteur ref ne sont pas        #
# répétées. L'usage entendu est que ref a la même valeur tout au long           #
# d'un certain bloc. Cela évite de comparer les valeurs de x à travers          #
# les blocs. Le paramètre acc permet de ne pas compter pour répétition les      #
# essais où il y a eu une erreur à l'essai -n. Pour cela, il faut que           #
# acc valle 0 dans le cas d'une erreur et 1 dans le cas d'un essai              #
# correct.                                                                      #
#===============================================================================#

makeRep <- function (x, ref, n, acc){
  v <- rep (NA, n)
  stimr <- c (v, x)
  stimp <- c (x, v)
  rep <- stimr==stimp
  rep <- rep[1:(length(rep)-n)]
  refr <- c (v, ref)
  refp <- c (ref, v)
  ref <- refr!=refp
  if (!missing (acc)){
    accr <- c (v, acc)
    accr <- accr[1:(length(acc)-n)]
    acc[accr==0] <- 1
    acc[accr==1] <- 0
    rep[acc==1] <- F
  }
  ref <- ref[1:(length(ref)-n)]
  rep[ref] <- F
  rep [is.na (rep)] <- F
  return (rep)
}


change.point <- function (x){
  N <- length (x)
  j <- 1:(N-1)
  Wj <- cumsum (rank (x[j]))
  K <- abs (2*Wj-(N+1)*j)
  Kmn <- max (K)
  m <- j[K==Kmn]
  n <- N-m
  o <- m*(N+1)/2
  Wm <- W[m]
  if (Wm > o) h <- -.5 else h <- .5
  z <- (Wm+h-o)/sqrt (n*o/6)
  structure (list (K=Kmn, m=m,z=z))
}


change.point.r <- function (x, bs=5000){
  N <- length (x)
  j <- 1:(N-1)
  Wj <- cumsum (rank (x[j]))
  K <- abs (2*Wj-(N+1)*j)
  Kmn <- max (K)
  m <- (j[K==Kmn])[1]
  count <- 0
  for (i in 1:bs){
    y <- sample (x)
    Wjy <- cumsum (rank (y[j]))
    Ky <- abs (2*Wjy-(N+1)*j)
    Kmny <- max (Ky)
    my <- (j[Ky==Kmny])[1]
    if ((Kmny>=Kmn)&(my<=m)) count <- count+1
  }
  p <- count/bs
  structure (list (K=Kmn, m=m,p=p))
}



# Construit un graphique d'interaction à deux facteurs, avec ou sans
# barres d'erreurs (calculées à la Stan), quand il s'agit d'un design
# avec mesures répétées. Le point crucial est ce qu'il faut lui
# fournir: la variable qu'on veut plotter et les facteurs, DONT LE
# FACTEUR ALEATOIRE. Ainsi, listoffactors doit être une liste dont le
# premier vecteur est typiquement "suj" et dont les deux suivants sont
# les facteurs dont on veut visualiser l'interraction. Ce qu'il faut
# noter, c'est qu'une fois qu'on a calculé pour chaque cellule du
# design et chaque "sujet" la valeur de la variable grâce à FUN, on
# aggrège cela avec une moyenne, codée en dur.
# Usage: inter.plot (rt, list (suj=suj, rule=rule, stim=stim))

# Donc, si je comprends bien, si tu donnes déjà une seule valeur par
# case, FUN ne sert à rien

inter.plot <- function (x, listoffactors,
                          FUN=median, errbars=T, dep=.05,
                          xlab,
                          ylab, main, at, ylim, labels=c(),...){
  c <- aggregate(x, listoffactors, FUN)
  if (missing (xlab)){xlab <- attributes (listoffactors)$names[2]}
  if (missing (ylab)){ylab <- ""}
  if (missing (main)){main <- ""}
  xpos <- levels (c[,2])
  if (missing (at)) {
    if (any (xpos==T)){at <- 0:1}
    else {at <- sort(numerize (xpos))}
  }
  e <- tapply(c$x, list (c[,2], c[,3]), mean)
  n <- dim (e)[2]
  if (n>2){leg <- c (1:2, 4:(n+1))} else {leg <- 1:2}
if (errbars==T){
  rtc <- rm.main.effect (c$x, c$suj, FUN)
  secells <- tapply(rtc,  list (c[,2], c[,3]), se)
  z <- dmatplot(at, e, secells, dep=dep, xlab=xlab, ylab=ylab,
           pch=leg, lty=leg, main=main, ylim=ylim, ...)}
else {
  z <- dmatplot(at, e, dep=dep, xlab=xlab, ylab=ylab,
           pch=leg, lty=leg, main=main, ylim=ylim, ...)}
               axis (2);axis (1,at=at, labels=labels)
structure (list (leg=leg, ylim=z$ylim, e=e, c=c))
}


 insert.value<-function(vec,newval,pos) {
  if(pos == 1) return(c(newval,vec))
  lvec<-length(vec)
  if(pos > lvec) return(c(vec,newval))
  return(c(vec[1:pos-1],newval,vec[pos:lvec]))
 }
 
permute<-function(elem) {
  if(!missing(elem)) {
   if(length(elem) == 2) return(matrix(c(elem,elem[2],elem[1]),nrow=2))
   last.matrix<-permute(elem[-1])
   dim.last<-dim(last.matrix)
   new.matrix<-matrix(0,nrow=dim.last[1]*(dim.last[2]+1),ncol=dim.last[2]+1)
   for(row in 1:(dim.last[1])) {
    for(col in 1:(dim.last[2]+1))
     new.matrix[row+(col-1)*dim.last[1],]<-insert.value(last.matrix[row,],elem[1],col)
   }
   return(new.matrix)
  }
  else cat("Usage: permute(elem)\n\twhere elem is a vector\n")
 } 

########################################
# EZ diffusion model (Wagenmakers 2005)#
########################################


logit <- function (p)
  {
    lp <- log(p/(1-p))
    return(lp)
  }


get.vaTer <- function (MRT, Pc, VRT,  s=.1)
  {
    s2 <- s^2
   # the default value for the scaling parameter s equals .1
    if (Pc == 0){
      cat("Oops, Pc == 0!\n")
      Pc <- .001}
    if (Pc == 0.5)
      cat("Oops, Pc == .5!\n")
    if (Pc == 1){
      cat("Oops, Pc == 1!\n")
      Pc <- .999
    }
   # If Pc equals 0, .5, or 1, the method will not work, and
# an edge-correction is required.
    L = logit(Pc)
    x = L*(L*Pc^2 - L*Pc + Pc - 0.5)/VRT
    v = sign(Pc-0.5)*s*x^(1/4)
                                        # this gives drift-rate
    a = s2*logit(Pc)/v
                                        # this gives boundary separation
    y = -v*a/s2
    MDT = (a/(2*v))*(1-exp(y))/(1+exp(y))
    Ter = MRT-MDT
                                        # this gives nondecision time
    return(list(v=v, a=a, Ter=Ter))
}


modalvalue <- function(x, na.rm=FALSE)
{
    x = unlist(x);
    if(na.rm) x = x[!is.na(x)]
    u = unique(x);
    n = length(u);
    frequencies = rep(0, n);
    for(i in 1:n)
    {
        if(is.na(u[i]))
        {
            frequencies[i] = sum(is.na(x))
        } else
        {
            frequencies[i] = sum(x==u[i], na.rm=TRUE)
        }
    }
    u[which.max(frequencies)]
}


mget.vaTer <- function(rt=NULL, acc=NULL, lof=NULL){
  nfact <- dim(data.frame(lof))[2]
  nam <- names(lof)
  fnam <- nam[1]
  if (nfact==1){
    d.lof <- list(data.frame(lof)[acc==1,])
    names(d.lof) <- fnam
  } else {
    d.lof <- data.frame(lof)[acc==1,]
  }
  x <- aggregate(rt[acc==1], d.lof, mean)
  names(x)[length(names(x))] <- "MRT"
  x <- cbind(x, Pc=aggregate(acc, lof, mean)$x)
  x <- cbind(x, VRT=aggregate(rt[acc==1], d.lof, var)$x)
  y <- c()
  for (i in 1:length(x[,1])){
    y <- rbind(y,unlist(get.vaTer(x$MRT[i], x$Pc[i], x$VRT[i])))
  }
  x <- cbind(x,y)
}

mget.vaTer4 <- function(rt=NULL, acc=NULL, lof=NULL){
  nfact <- dim(data.frame(lof))[2]
  nam <- names(lof)
  fnam <- nam[1]
  if (nfact==1){
    d.lof <- list(data.frame(lof)[acc==1,])
    names(d.lof) <- fnam
  } else {
    d.lof <- data.frame(lof)[acc==1,]
  }
  x <- aggregate(rt[acc==1], d.lof, mean)
  names(x)[length(names(x))] <- "MRT"
  x <- cbind(x, Pc=aggregate(acc, lof, mean)$x)
  x <- cbind(x, VRT=aggregate(rt[acc==1], d.lof, var)$x)
  x$Pc = (x$Pc-.25)*50/75+.5001
  y <- c()
  for (i in 1:length(x[,1])){
    y <- rbind(y,unlist(get.vaTer(x$MRT[i], x$Pc[i], x$VRT[i])))
  }
  x <- cbind(x,y)
}


## 2013 - 02 - 25


get.sigmaMuTau <- function(rt, lof){
  
## Cette fonction applique la fonction "timefit" du package retimes
## sur les temps de réponses dans chaque cellule de lof. Elle retourne
## un data.frame du type retourné par aggregate, avec en colonne, la
## médiane des temps de réponses, et les paramètres de l'ex-gaussienne
## ajustée sur lesdits temps de réponses.

  df <- data.frame(rt=rt, lof)
  KK <- evalq(by(rt, lof, simplify=FALSE, FUN=function(x){as.vector(timefit(x)@par)}), df)
  KM <- matrix(unlist(KK), byrow=TRUE, ncol=3)
  e <- evalq(aggregate(rt, lof, median), df)
  dimnames(e)[[2]][length(dimnames(e)[[2]])] <- "median"
  dimnames(KM)[[2]] <- c("mu", "sigma", "tau")
  return(cbind(e, KM))
}




###############

## Pour fitter des sigmoïdes de fonctions psychométriques:



fitbin <- function(dur, est, n=10, plot=TRUE, xlim=c(200,1100), ylim=c(0,1), main='', add=FALSE, method='quantiles', cex=.6,
                   lwd=2){
  if (method=='quantiles'){
    l <- bin.count(dur, est, n)
  }
  else {
    l <- bin.count.fixed(dur, est, n)
  }
  resp.mat <- matrix(c(l$t1, l$t0), ncol=2)
  phat <- as.vector(l$t1)/(as.vector(l$t0)+as.vector(l$t1))
  cnt <- l$x;z <- seq(min(dur), max(dur), by=1)
  l <- glm(resp.mat~cnt, family=binomial)
  if (plot){
    if (add){
    points(cnt, phat, col=3, cex=cex, pch=2)
    lines(z, predict(l, data.frame(cnt=z), type='response'), col=3, lty=2, lwd=lwd)
    lines(c(0,10000), c(.5, .5), lty=2, lwd=.8)
  }
    else
      {
        plot(cnt, phat, xlim=xlim, ylim=ylim, main=main, col=2, cex=cex, pch=2)
        lines(z, predict(l, data.frame(cnt=z), type='response'), col=2, lty=2, lwd=lwd)
        lines(c(0,10000), c(.5, .5), lty=2, lwd=.8)
      }
  }
  list(pse=dose.p(l), coef=l$coefficients, p.values=summary(l)$coefficients[,4], cnt=cnt, phat=phat)
}


     pchShow <-
       function(extras = c("*",".", "o","O","0","+","-","|","%","#"),
                cex = 3, ## good for both .Device=="postscript" and "x11"
                col = "red3", bg = "gold", coltext = "brown", cextext = 1.2,
                main = paste("plot symbols :  points (...  pch = *, cex =", cex,")"))
       {
         nex <- length(extras)
         np  <- 26 + nex
         ipch <- 0:(np-1)
         k <- floor(sqrt(np))
         dd <- c(-1,1)/2
         rx <- dd + range(ix <- ipch %/% k)
         ry <- dd + range(iy <- 3 + (k-1)- ipch %% k)
         pch <- as.list(ipch) # list with integers & strings
         if(nex > 0) pch[26+ 1:nex] <- as.list(extras)
         plot(rx, ry, type="n", axes = FALSE, xlab = "", ylab = "", main = main)
         abline(v = ix, h = iy, col = "lightgray", lty = "dotted")
         for(i in 1:np) {
           pc <- pch[[i]]
           ## 'col' symbols with a 'bg'-colored interior (where available) :
           points(ix[i], iy[i], pch = pc, col = col, bg = bg, cex = cex)
           if(cextext > 0)
               text(ix[i] - .3, iy[i], pc, col = coltext, cex = cextext)
         }
       }


# Fonctions de Vincent de Gardlle pour les anovas

######
# pour catégoriser et labeliser les valeurs d'un vecteur avec des seuils, ex des pvalues
######
pvseuils <- c(Inf, 0.1,0.05,0.01,0.001)
pvlabels <- c(" ",".","*","**","***")
dostars <- function(vect,seuils=pvseuils,labels=pvlabels) {
newvect <- vect
for (i in 1:length(seuils)) {newvect[vect<seuils[i]] <- labels[i]}
dostars <- newvect
}

######
# graph adhoc (modifiable) pour plots avec 3 facteurs
######
mygraph2 <- function(measure,fac1,fac2,absc,mycols,ylims,leg="TRUE") {
# for 3-factors plots (absc is the factor used for the abscissa)
# measure, fac1, fac2, absc, should be same size, for instance taken from the same dataframe
# ylims = c(lower limit, upper limit)
# mycols = color vector, e.g. mycols <- c("dark green","green","dark blue","blue")
icol <- c(0)
mylegend <- mycols
for (ifac1 in levels(fac1)) {
for (ifac2 in levels(fac2)) {
X <- fac1==ifac1 & fac2==ifac2
icol <- icol + 1
graph<- tapply(measure[X],list(absc[X]),mean)
plot(graph, type = "b", col = mycols[icol], ylim=ylims)
par(new =TRUE)
mylegend[icol] <- paste(ifac1,ifac2)
}}
par(new=FALSE)
if (leg) {legend("topright", fill = mycols, legend = mylegend,inset=0.02)}
print(rbind(mylegend,mycols))
return(mylegend)
}

######
# pour sortir les resultats de l'anova au format standard,
# avec F val (et df), p val, et eta sq.
######

    ######
    # pour sortir les resultats de l'anova au format standard
    myprintaov <- function(a) {
    # a est un objet anova, ex: a <- aov(rt~condition1 + Error(suj/condition1),data=priming)
    err_lay <- names(summary(a))

    for (ierr in err_lay) {
    cur_lay <- eval(parse(text=paste("summary(a)$\"",ierr,"\"",sep="")))[[1]]

    nb_fact_lay <- nrow(cur_lay)-1
    nm_fact_lay <- rownames(cur_lay)[1:nb_fact_lay]

    df_cur_lay  <-  cur_lay$Df[nrow(cur_lay)]
    SS_cur_lay <- sum(cur_lay$"Sum Sq")

    for (ifact in 1:nb_fact_lay) {
    df_cur_fact <- cur_lay$Df[ifact]
    SS_cur_fact <- cur_lay$"Sum Sq"[ifact]
    Fv_cur_fact <- round(cur_lay$"F value"[ifact],digits=2)
    pv_cur_fact <- cur_lay$"Pr(>F)"[ifact]
    pvseuils <- c(Inf, 0.1,0.05,0.01,0.001)
    pv_thr  <- last(pvseuils[pv_cur_fact < pvseuils])
    etasquare <- round(SS_cur_fact/SS_cur_lay,digits=3)
    print(paste(nm_fact_lay[ifact]," (F(", df_cur_fact, ",", df_cur_lay, ")=", Fv_cur_fact, ", p<", pv_thr, ", eta-sq=", etasquare, ")", sep=""))
    }}}



     pchShow <-
       function(extras = c("*",".", "o","O","0","+","-","|","%","#"),
                cex = 3, ## good for both .Device=="postscript" and "x11"
                col = "red3", bg = "gold", coltext = "brown", cextext = 1.2,
                main = paste("plot symbols :  points (...  pch = *, cex =",
                             cex,")"))
       {
         nex <- length(extras)
         np  <- 26 + nex
         ipch <- 0:(np-1)
         k <- floor(sqrt(np))
         dd <- c(-1,1)/2
         rx <- dd + range(ix <- ipch %/% k)
         ry <- dd + range(iy <- 3 + (k-1)- ipch %% k)
         pch <- as.list(ipch) # list with integers & strings
         if(nex > 0) pch[26+ 1:nex] <- as.list(extras)
         plot(rx, ry, type="n", axes = FALSE, xlab = "", ylab = "",
              main = main)
         abline(v = ix, h = iy, col = "lightgray", lty = "dotted")
         for(i in 1:np) {
           pc <- pch[[i]]
           ## 'col' symbols with a 'bg'-colored interior (where available) :
           points(ix[i], iy[i], pch = pc, col = col, bg = bg, cex = cex)
           if(cextext > 0)
               text(ix[i] - 0.3, iy[i], pc, col = coltext, cex = cextext)
         }
       }


################################################################################################
## mAFC dprime calculations. Ceci est basé sur JE Keith Smith, 1982, lui-même basé sur Hacker ##
## and Ratcliff 1979, et cité, ainsi que le précédent, par MacMillan and Creelman. Bref, rien ##
## que du beau monde. Bizarre que cette formule ne soit pas plus répandue...                  ##
##                                                                                            ##
## NOTE CRITIQUE: CELA SUPPOSE PAS DE BIAIS!!                                                 ##
################################################################################################

 

## gen.dprime2 <- function (m, pc){
##   l <- log(m-1)
##   A <- (-4 + sqrt(16+25*l))/3
##   B <- sqrt((l+2)/(l+1))
##   dprime <- A - B*qnorm(pc, lower.tail=F)
##   return(dprime=dprime)
## }

gen.dprime <- function (m, pc){
  km <- 0.86 - 0.085 * log(m-1)
  dprime <- km * log((m-1)*pc/(1-pc))
  return(dprime=dprime)
}





# FONCTIONS POUR LES ROCs EMPIRIQUES

# Les fonctions suivantes retournent chacune un type différent de ROC,
# tracent la courbe et renvoient l'aire sous la courbe et la première
# diagonale (auc) ainsi qu'une estimation du biais (B).

# Il y avait plusieurs erreurs dans le code dans les versions avant
# 2014 - 09 - 07.

# roc1 fait cela pour la ROC de premier ordre en utilisant la
# confiance (trustR) pour fixer a posteriori des critères de réponses
# différents (cf MacMillan et Creelman, chap.3). A priori, on n'en n'a
# pas besoin, mais c'est là pour mémoire.

# roc2 fait ça pour la ROC de second ordre, aussi en utilisant la
# confiance. Ça permet d'avoir une estimation empirique de la
# sensibilité de second ordre non biaisée. C'est ce que font Fleming
# et al. 2010 ainsi que Kolb et Braun 1995.

# rocrt fait ça pour les temps de réponses estimés (je crois que c'est
# tout à fait nouveau comme type d'analyse.) La logique est la
# suivante: on partage les temps de réponses (objectifs et estimés) en
# deux quantiles (par exemples les 10% plus rapides vs les 90 % les
# plus lents; ou bien les 90% plus rapides et les 10% les plus lents).
# Pour chacun de ces critères, un hit est un essai tel que les temps
# de réponses estimés et objectifs soient tous les deux rapides. Le
# pourcentage de hits est obtenu en divisant par le nombre d'essais
# objectivement rapides (selon ce critère). Une fausse alarme est un
# essai estimé rapide alors qu'il est objectivement lent; le
# pourcentage de fausses alarmes est obtenu en divisant par le nombre
# d'essais objectivement lents selon ce critère. On bouge le critère
# et on obtient ainsi plusieurs couples (hits, fas), qui nous
# permettent de tracer une ROC et une auc. C'est une estimation de la
# sensibilité métacognitive au temps de réponse, non biaisée et
# surtout non paramétrique puisque purement ordinale (ce que n'est pas
# la corrélation).

# Voir des exemples d'utilisation dans behav.r de l'analyse des
# données de la deuxième manip avec mesure du cortisol de Gabriel.

roc2 <- function(conf, acc, n=7, plot=TRUE, add=FALSE, col=NULL, title='', cex=2, lwd=1){
  stopifnot(require(pracma))
  s <- seq(1, 0, length.out=n)
  # A second order hit is: correct with high confidence. A second
  # order false alarm is: error with high confidence. To get the
  # rates, you divide by the overall number of correct and incorrect
  # trials, respectively. Since you want the ROC, you do that for a
  # number of confidence levels. Notice that conservative points (near
  # (0,0) in the ROC square) come from high confidence criteria (near
  # 1 on the confidence scale.) This is a bit counter intuitive, and
  # is the reason behind ordering confidence criteria from 1 to 0.
  fhit <- function(x) {sum(acc==1 & conf >= x)/sum(acc==1)}
  ffas <- function(x) {sum(acc==0 & conf >= x)/sum(acc==0)}
  hits <- c(0, unlist(lapply(s, fhit)))
  fas <- c(0, unlist(lapply(s, ffas)))
  # We use the geometrical utilities in package "pracma"; truncSquare
  # is a large polygon with the minor diagonal of the ROC space, that
  # we use to get the intersection.
  truncSquare <- matrix(c(-1, -1, 0, 1, 1, -1, 1, 1, 0, -1), byrow=TRUE, nrow=2)
  AUC <- rbind(fas, hits)
  biasXY <- poly_crossings(AUC, truncSquare)[1,]
  cutPoint <- AUC[1,]<biasXY[1]
  AUCL <- cbind(rbind(fas[cutPoint], hits[cutPoint]), cbind(biasXY, .5))
  AUCR <- cbind( cbind(.5, biasXY),rbind(fas[!cutPoint], hits[!cutPoint]))
  # We reverse the order of the points to get positive areas, because
  # pracma uses the trigonometrical convention.
  auc <- polyarea(rev(AUC[1,]), rev(AUC[2,]))
  aucL <- polyarea(rev(AUCL[1,]), rev(AUCL[2,]))
  aucR <- polyarea(rev(AUCR[1,]), rev(AUCR[2,]))
  # B is the bias: 0 means unbiased (as much area above ("right") and
  # below ("left") of the minor diagonal); B positive means more area
  # above ("right", near the (1,1) point in the ROC space, that is
  # many hits and many FAs---ie LIBERAL); B negative means more area
  # below ("left": CONSERVATIVE). Since it's a ratio we report the
  # log.
  B <- log(aucR/aucL)
  if (plot){
      if (!add){
          main=paste(title, 'AUC: ', round(auc, 3), 'B: ', round(B, 3))
          plot(fas, hits, type='b', cex=cex, lty=1, xlim=c(0, 1), ylim=c(0, 1),
               xlab='', ylab='', col=col,
               main=main, pch=19)
          lines(0:1, 0:1, lty=2, lwd=lwd)
      } else {
          lines(fas, hits, type='b', cex=cex, lty=1, col=col, pch=19)
      }
  }
  return(list(auc2=auc, B2=B, fas=fas, hits=hits))
}

roc1 <- function(signedConf, stim, n=7, plot=TRUE, add=FALSE, col=NULL){
  stopifnot(require(pracma))
  s <-  seq(-1, 1, length.out=n)
  # signedConf is a fullscale confidence scale, from -1 to 1. Meaning
  # that it combines the response and confidence (response for one
  # stimulus goes from -1 to 0 and for the other from 0 to 1). Put
  # your first confidence criterion at 0: if stim == -1 and signedConf
  # < 0, you have a hit; but if stim = -1, it's a false alarm. Now,
  # generalize for other criteria...
  fhit <- function(x) {sum(stim==-1  & signedConf <= x )/sum(stim==-1)}
  ffas <- function(x) {sum(stim== 1  & signedConf <= x )/sum(stim==1)}
  hits <- c(0, unlist(lapply(s,fhit)) )
  fas <- c(0, unlist(lapply(s,ffas)) )

  # Geometry is done the same way as in roc2
  truncSquare <- matrix(c(-1, -1, 0, 1, 1, -1, 1, 1, 0, -1), byrow=TRUE, nrow=2)
  AUC <- rbind(fas, hits)
  biasXY <- poly_crossings(AUC, truncSquare)[1,]
  cutPoint <- AUC[1,]<biasXY[1]
  AUCL <- cbind(rbind(fas[cutPoint], hits[cutPoint]), cbind(biasXY, .5))
  AUCR <- cbind( cbind(.5, biasXY),rbind(fas[!cutPoint], hits[!cutPoint]))
  auc <- polyarea(rev(AUC[1,]), rev(AUC[2,]))
  aucL <- polyarea(rev(AUCL[1,]), rev(AUCL[2,]))
  aucR <- polyarea(rev(AUCR[1,]), rev(AUCR[2,]))
  
  B <- log(aucR/aucL)
  if (plot){
    if (!add){
      plot(fas, hits, type='b', lty=1, xlim=c(0, 1), ylim=c(0, 1), col=col,
           main=round(auc, 3))
    } else {
      lines(fas, hits, type='b',  lty=1, col=col)
    }
  }
  return(list(auc1=auc, B1=B))
}

rocrt <- function(rt, est, n=6, plot=TRUE, add=FALSE, col=1){
  stopifnot(n%%2==0); mid <- n/2+1
  N <- length(rt)
  qRtN <- cut2(rt, g=n, levels.mean=TRUE)
  qEstN <- cut2(est, g=n, levels.mean=TRUE)

  qRt <- as.numeric(qRtN)
  qEst <- as.numeric(qEstN)
  hits <- c(); fas <- c()
  for (i in 1:n){
    hits <- c(hits, sum(qEst<=i & qRt<=i)/sum(qRt<=i))
  }
  for (i in 1:(n-1)){
    fas <- c(fas, sum(qEst<=i & qRt>i)/sum(qRt>i))
  }
  fas <- c(fas, 1)
  hits <- c(0, hits); fas <- c(0, fas)
  fasp1 <- c(NA, fas); hitsp1 <- c(NA, hits); fas <- c(fas, NA); hits <- c(hits, NA)
  auc <- sum((hits-fasp1)^2-(hitsp1-fas)^2, na.rm=TRUE)/4
  aucL <- sum((hits[1:mid]-fasp1[1:mid])^2-(hitsp1[1:mid]-fas[1:mid])^2, na.rm=TRUE)/4
  aucR <- sum((hits[(mid+1):(n+1)]-fasp1[(mid+1):(n+1)])^2-(hitsp1[(mid+1):(n+1)]-fas[(mid+1):(n+1)])^2, na.rm=TRUE)/4
  B <- log(aucR/aucL)
  if (plot){
    if (!add){
      plot(fas, hits, type='b', cex=.5, lty=1, xlim=c(0, 1), ylim=c(0, 1),
           main=round(auc, 3))
      lines(0:1, 0:1, lty=2)
    } else {
      lines(fas, hits, type='b', cex=.5, lty=1, col=col)
    }
  }
 return(list(auc=auc, B=B))
}



drawCh <- function(M){
  # Draws the convex hull around a 2XN matrix of points
  ch <- chull(M)
  lines(M[ch,])
  l <- tail(ch, 1)
  lines (c(M[l, 1], M[ch[1],1]), c(M[l,2], M[ch[1],2]))
}

peeling <- function(M, alpha=.05){
  
  # Draws the convex hull around a 2XN matrix of points, after
  # removing alpha percent of those. It aims at being as close to this
  # percentage as possible.

  total <- length(M[,1])
  stopifnot( total>3 )
  removed <- 0
  ch <- c()
  while(removed / total < alpha ) {
    prevCh <- ch; prevM <- M
    prevRemoved <- removed
    ch <- chull(M)
    removed <- removed + length(ch)
    M <- M[-ch,]
  }
  if ( abs ((removed/total)-alpha) > abs ((prevRemoved/total)-alpha) ){
    M <- prevM
    drawCh(M)
    print(removed)
    removed <- prevRemoved
    if (prevRemoved == 0){stop("first peeling too large! increase alpha!")}
  }
  drawCh(M)
  return(list(ch=ch, error=removed/length(M[,1])))
}


readFastDM <- function(file, dv='v', factors, flevels){

    ##    lit la sortie d'un script fast-dm, et met le résultat dans
    ##    un data.frame. Voir un exemple d'application dans le
    ##    tutorial sur fast-dm.

    fdm <- read.table(file=file, header=TRUE)
    cols <-  dimnames(fdm)[[2]]
    assign(dv, unlist(fdm[,grep(paste('^',dv, '_', sep=''), cols)]))
    i <- 0
    reps <- prod(unlist(lapply(flevels, function(s){length(s)})))
    df <- data.frame(rep(fdm$dataset, reps))    

    for (f in factors){
        i <- i+1
        assign(f, rep(flevels[[i]][1], length(get(dv))))
        for (j in 2:length(flevels[[i]])){
            tmp <- get(f)
            tmp[grep(paste('_', flevels[[i]][j], sep=''), names(get(dv)))] <- flevels[[i]][j]
            assign(f, tmp)
        }
        df <- cbind(df, get(f))
    }
    df <- cbind(df, get(dv))
    names(df) <- c('suj', factors, dv);row.names(df) <- NULL
    df$suj <- factor(df$suj)
    return(df)
}

readFastDMComplete <-  function(file, dvs, factors, flevels){
    
    ## Si tu as utilisé exactement les mêmes facteurs et niveaux pour
    ## tous les paramètres de diffusion, tu peux utiliser cette
    ## fonction pour tout récupérer en même temps.

    L <- lapply(dvs, function(dv){readFastDM(file=file, dv=dv, factors=factors, flevels=flevels)})
    DF <- L[[1]]
    for (i in 2:length(L)){
        DF <- merge(DF, L[[i]])
    }
    return(DF)
}
    

readFastDMFixed <- function(file, dvs){

## Reads the parameters that haven't been fitted between conditions.
## Usually "sz", "st0" and "sv".

fdm <- read.table(file=file, header=TRUE)
df <- data.frame(suj=fdm$dataset)
for (dv in dvs){
    cols <-  dimnames(fdm)[[2]]    
    assign(dv, unlist(fdm[,grep(paste('^',dv, sep=''), cols)]))
    
    df <- cbind(df, get(dv))
}
names(df) <- c('suj', dvs);row.names(df) <- NULL
df$suj <- factor(df$suj)
return(df)
}

## The following function gets the mean of the quantiles along which
## you cut a variable. Thus:

## l  <- cut(data$var, breaks=quantile(data$var, probs=seq(0, 1, length.out=4)))
## data$catVar <- getQuantiles(l)

getMeanOfQuantiles <- function(liste){
    L <- lapply(l, function(m){
        ll <- strsplit(gsub("\\(|\\]", '', as.character(m)), ',', fixed=TRUE)
        ll <- mean(as.numeric(unlist(ll)))
        return(ll)})
    return(unlist(L))
}

## 2017 - 05 - 18

joinMCMCLists <- function(L){
    ## We check that all elements have the same number of chains. This
    ## function is intended to join updates!! (that is, run on the
    ## same model, with the same seed, etc.)
    stopifnot(length(unique(lapply(L, length)))==1)
    nChains <- length(L[[1]])
    thinnings <- unique(lapply(L, thin))
    stopifnot(length(thinnings)==1)
    th <- unique(unlist(thinnings))
    st <- min(unlist(lapply(L, start)))
    en <- max(unlist(lapply(L, end)))
    LL <- mcmc.list(lapply(1:nChains, function(i){
                               M <- mcmc(do.call(rbind, lapply(L, function(l){l[[i]]})),
                                         thin=th, start=st, end=en)
                    }))
}


etas <- function(l){
    summ <- summary(l)
    for (i in 2:(length(summ))){
        print(names(summ)[i])
        eta <- summ[[i]][[1]]$Sum[1]/(summ[[i]][[1]]$Sum[1] + summ[[i]][[1]]$Sum[2])
        print(eta)
      }
}



print('finished sourcing R_Lib_Sackur.r!!')

