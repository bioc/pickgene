#####################################################################
##
## $Id: newton.R,v 1.1 2003/11/25 14:50:35 jgentry Exp $
##
##     Copyright (C) 1999, 2000 Michael A. Newton.
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## The text of the GNU General Public License, version 2, is available
## as http://www.gnu.org/copyleft or by writing to the Free Software
## Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
##
###############################################################################
##
##Code in ftp.biostat.wisc.edu/pub/newton/Arrays/code/
## was used to implement the calculations reported
##in Newton et al. 1999. On differential variability of expression ratios:
##Improving statistical inference about gene expression changes from 
##microarray data.  Submitted to J. Comp. Biol., 11/1999.
##
##See www.stat.wisc.edu/~newton/ for further information.
##
##The files all use data read in by code in the `read' directory.
##At the moment, the raw data files are unavailable. To implement
##calculations, simply make sure that intensity measurements get read
##into two vectors, ``xx'' and ``yy'' of length equal to the number of
##spots on the microarray.
##
##Briefly, the files in this directory do the following:
##
##oddsplot        plots the odds of change, Fig. 4
##                  (uses fits from em.ggb)
##
##em.ggb              fits the Gamma-Gamma=Bernoulli model via EM
##
##pmarg           computes the profile loglikelihood for p
##
##fitgg           fits the Gamma-Gamma model to one array
##
##rankgene            compares the ranking of genes by the naive procedure
##                  and the empirical bayes procedure; uses fits stored
##                  in ../results/fits.gg (and makes Fig. 2)
##
##datplot         creates scatterplots of intensity measurements
##
##shrinkplot      plots Fig. 1, showing shrinkage
##
##s.check0          compares marginal histograms to fitted margins (Fig. 5)
##
##s.check1          diagnostic check (Fig. 6)
##
##s.check2          looks at some of the changed genes
##
#######################################################################
do.oddsplot <- function(data,
                        main = substitute( data ),
                        theta = c(2,2,2,.4),
                        col = NULL,
                        xlab = conditions[1], ylab = conditions[2],
                        redo = missing( theta ),
                        conditions = c("Cy3","Cy5"),
                        identifier = "identifier", ... ) {
  if( redo )
    theta <- em.ggb(data[[conditions[1]]], data[[conditions[2]]],
                    theta, theta[1:3], print = TRUE )
  lod <- oddsplot( data[[conditions[1]]], data[[conditions[2]]], theta,
                  xlab = xlab, ylab = ylab,
                  main = main, col = col, ... )
  if( ncol( data ) > 2 )
    probes <- data[[identifier]]
  else
    probes <- seq( nrow( data ))
  probes <- lodprobes( data[[conditions[1]]], data[[conditions[2]]], theta,
                      lod, probes, col = col )
  print( probes )
  invisible( list( theta = theta, lod = lod, probes = probes ))
}
#######################################################################
normal.richmond <- function( foo = read.table( "../data/mn2.csv", header = TRUE,
                               sep = "," ),
                            channel = "BP109CH" )
{
  ## Normalize to average intensity first using Richmond et al.

  nspot <- nrow(foo)

  ## Background adjustment using channel (very simple)
  x <- foo[[ paste( channel, "1I", sep = "" ) ]] -
    foo[[ paste( channel, "1B", sep = "" ) ]]
  y <- foo[[ paste( channel, "2I", sep = "" ) ]] -
    foo[[ paste( channel, "2B", sep = "" ) ]]

  ## Normalization
  ## Rescale to help with underflow problem 10^5 (does not affect shape params)
  x <- 100000 * x / sum( x[x>0] )
  y <- 100000 * y / sum( y[y>0] )

  ok <- x>0 & y>0
  list(xx = x[ok], yy = y[ok] )
}
#######################################################################
chen.poly <- function(cv,err=.01)
{
  ## part of table 2 from Chen et al
  bar <- rbind( c(.979, -2.706, 2.911, -2.805 ),
               c(.989, 3.082, -2.83, 28.64),
               c(.9968, -3.496,4.462, -5.002),
               c( .9648,4.810,-15.161,78.349) )
  if( err==.05 ) {
    coef <- bar[1,]
    tmp1 <- cv^3*coef[4] + cv^2*coef[3]+cv*coef[2] + coef[1]
    coef <- bar[2,]
    tmp2 <- cv^3*coef[4] + cv^2*coef[3]+cv*coef[2] + coef[1]
  }
  if( err==.01 ) {
    coef <- bar[3,]
    tmp1 <- cv^3*coef[4] + cv^2*coef[3]+cv*coef[2] + coef[1]
    coef <- bar[4,]
    tmp2 <- cv^3*coef[4] + cv^2*coef[3]+cv*coef[2] + coef[1]
  }
  c(tmp1,tmp2)
}
#######################################################################
fitgg <- function( xx, yy, start = c(10,1,1) )
## green in xx, red in yy
{
  ## Fits the Gamma-Gamma model 

  bar <- nlminb( start=start, objective=nloglik, lower=c(1,0,0),
		xx = xx, yy = yy )
  fits <- c( bar$par, length(xx) )
  names( fits ) <- c("aa","a0","nu","n")
  fits
}
#######################################################################
# following uses R's nlm() as surrogate for Splus's nlminb()
# comment out or delete when used in Splus
# See notes on nloglik and nploglik below when used in R
#######################################################################
nlminb <- function( start=c(10,1,1), objective, lower=c(1,0,0), xx, yy, zz,
                   use.optim = FALSE )
{
  ## kludge to make xx, yy and zz global to nloglik or nploglik
  if( !missing( xx ))
    assign( ".fit.xx", xx, pos = 1 )
  if( !missing( yy ))
    assign( ".fit.yy", yy, pos = 1 )
  if( !missing( zz ))
    assign( ".fit.zz", zz, pos = 1 )

  ## also has optim which does
  if( use.optim )
    theta <- optim( start, objective, lower = lower,
                   method = "L-BFGS-B" )$par
  else {    
    ## R has routine nlm, which does not take care of bounds
    ## so we just redefine parameters
    if( !missing( lower )) {
      for( i in seq( length( start ))) {
        if( lower[i] == 0 )
          start[i] <- log( start[i] )
        if( lower[i] == 1 )
          start[i] <- log( log( start[i] ))
      }
    }
    theta <- nlm( objective, start )$estimate

    ## and now we backtransform the parameters
    if( !missing( lower )) {
      for( i in seq( length( start ))) {
        if( lower[i] == 0 )
          theta[i] <- exp( theta[i] )
        if( lower[i] == 1 )
          theta[i] <- exp( exp( theta[i] ))
      }
    }
  }
  theta
}
#######################################################################
lod.ggb <- function(x,y,theta)
{
  ## Log_(10) posterior odds
  ## x = channel 1 intensity
  ## y = channel 2 intensity

  ## theta = (aa,a0,nu,pp)
  aa <- theta[1]; a0 <- theta[2]
  z0 <- y0 <- x0 <- theta[3]
  pp <- theta[4]
  tmp <- log( pp ) - log(1-pp) +
    a0*( log(x0) + log(y0) - log(z0) ) +
      (2*aa+a0)*log(x+y+z0) -
        (aa+a0)*( log(x+x0) + log(y+y0) ) +
          2*lgamma(aa+a0) - lgamma(a0) - lgamma(2*aa+a0)
  tmp / 2.3
}
#######################################################################
loglik <- function(theta,xx,yy)
{
  ## Returns loglikelihood for observed data
  ## xx,yy are intensities in the two channels

  ## theta=(aa,a0,nu,p)
  aa <- theta[1]; a0<-theta[2]; nu<-theta[3]

  n <- length(xx)

  ## p_0(r,g) (with common factor (rg)^(a-1) removed

  lp0 <- lgamma(2*aa+a0) + a0*log(nu) - 2*lgamma(aa) - lgamma(a0) -
    (2*aa+a0)*log( xx+yy+nu )

  ## p_a(r,g)
  lpa <-  2*(lgamma(aa+a0)-lgamma(aa)-lgamma(a0)) +
    + 2*a0*log(nu) - (aa+a0)*log( (xx+nu)*(yy+nu) ) 

  ll <- (aa-1)*log(xx*yy) + log( theta[4]*exp(lpa) +
                                (1-theta[4])*exp(lp0) )
  return(sum(ll))
}
#######################################################################
nloglik <- function( theta, xx = .fit.xx, yy = .fit.yy )
{
  ## theta=(log(log(aa)),log(a0),log(nu))
  ## uncomment the following two lines if used in R
  theta <- exp( theta )
  theta[1] <- exp( theta[1] )
  aa <- theta[1]; a0<-theta[2]; x0<-theta[3]; y0<- theta[3]

  n <- length(xx)

  ll <- 2*n * ( lgamma(aa+a0) - lgamma(aa) - lgamma(a0) )
  ll <- ll + n*a0 * ( log(x0)+log(y0) ) + (aa-1) * sum( log(xx)+log(yy) )
  (aa+a0) * sum( log(x0+xx) + log(y0+yy) ) - ll
}
#######################################################################
nploglik <- function( theta, xx= .fit.xx, yy = .fit.yy, zz = .fit.zz )
{
  ## xx,yy are intensities in the two channels; zz=P(b!=c|xx,yy)
  ## (I'll separately optimize pp=P(zz=1); hence npl.. for partial loglik

  ## theta=(log(log(aa)),log(a0),log(nu))
  ## uncomment the following two lines if used in R
  theta <- exp( theta )
  theta[1] <- exp( theta[1] )

  aa <- theta[1]; a0<-theta[2]; x0<-theta[3]; y0<- theta[3];
  z0 <- theta[3]
  n <- length(xx)

  ## Complete data loglikelihood
  sumzz <- sum( zz )
  lgaa <- lgamma( aa )
  lga0 <- lgamma( a0 )
  ll <- (aa-1) * sum( log(xx) + log(yy) ) +
    sumzz * 2 * ( lgamma(aa+a0) - lgaa - lga0 ) +
      sumzz*a0*(log(x0)+log(y0)) +
        (n-sumzz) * ( lgamma(2*aa+a0) - 2 * lgaa - lga0 ) +
          (n-sumzz) * a0 * log(z0) -
            (aa+a0) * sum( zz * ( log(x0+xx) + log(y0+yy) ) ) -
              (2*aa+a0) * sum( (1-zz) * ( log(z0+xx+yy) ) )
  -ll
}
#######################################################################
rankgene <- function( xx, yy, fits = fitgg( xx, yy ))
{
  ## Look at effect on rank of the shrinkage
  ## Shrinkage factors from fits.gg
  xhat <- xx + fits[3]
  yhat <- yy + fits[3]
  eps <- runif( length(xx) ) * .00001   # randomize a bit to break ties
  
  r1a  <- rank( abs( log(xx/yy) + eps ) )  # raw ranking (most change, either way)
  r2a  <- rank( abs( log(xhat/yhat) + eps ) ) # Bayes ranking

  ## Look at top 100 genes most changed by raw ranking
  counta <- rep(NA,100)
  na <- length(r1a)
  for( i in 1:100 ) {
    ind <- (1:na)[ r1a <= i ]   # highly ranked by raw method
    counta[i] <- sum( r2a[ind] <= i )  # how many similarly ranked
  }
  list( count = counta, r1 = r1a, r2 = r2a )
}
#######################################################################
em.ggb <- function( x, y, theta = c(2,2,2,.4), start = c(2,1.2,2.7),
                 pprior = 2, printit = FALSE, tol = 1e-9, offset = 0 )
{
### Fit Gamma/Gamma/Bernoulli model (equal marginal distributions)
### 
### Model:
### spot intensities x ~ Gamma(a,b); y ~ Gamma(a,c)
### w.p. p,	b=c,		common value ~ Gamma(a0,nu)
###     w.p. 1-p,	b != c, 	values ~ Gamma(a0,nu)
###     all independent

  tmp <- x > -offset & y > -offset
  x <- x[tmp] + offset
  y <- y[tmp] + offset
  if( any( !tmp ))
    warning( paste( sum( !tmp ), "probes dropped with values below", offset ))
  rm( tmp )
  n <- length(x)
  if( pprior ) {
    ## kludge to make x and y global to nploglik
    assign( ".fit.xx", x, pos = 1 )
    assign( ".fit.yy", y, pos = 1 )
  }

  ## EM algorithm

  ## starting value
  notdone <- TRUE
  iter <- 1
  while( notdone ) {
    aa <- theta[1]; a0<-theta[2]; x0<-theta[3]; y0<- theta[3];
    z0 <- theta[3]; pp <- theta[4]

    ## E-step 
    tmp <- log( pp ) - log(1-pp) +
      a0*( log(x0) + log(y0) - log(z0) ) +
        (2*aa+a0)*log(x+y+z0) -
          (aa+a0)*( log(x+x0) + log(y+y0) ) +
            2*lgamma(aa+a0) - lgamma(a0) - lgamma(2*aa+a0)
    zz <- 1/( 1 + exp(-tmp) )

    ## M-step
    fit <- nlminb( start=start, objective=nploglik,
                  lower=c(1,0,0), zz=zz )

    ## check tolerance
    chk <- sum(( theta[1:3] - fit$parameter )^2 )
    
    ## Add a prior on pp
    theta[1:3] <- fit$parameter

    ## Beta hyperparameter for p
    if( pprior )
      theta[4] <- ( pprior + sum( zz ) ) / ( 2 * pprior + n )
    if( printit )
      print(round(theta,4) )
    iter <- iter + 1
    notdone <- (chk > tol) & (iter<100)
  } 
  theta
}
#######################################################################
pmarg <- function( xx, yy, theta = c(2.75,1.37,4.12), nsupp = 20 )
{
# This file gets a profile loglikelihood for the mixing rate p

  ## kludge to make xx and yy global to nploglik
  assign( ".fit.xx", xx, pos = 1 )
  assign( ".fit.yy", yy, pos = 1 )

  ## support for heat-shock example
  psupp <- seq( .0001, .2, length = nsupp )
  lprof <- array( NA, 5, nsupp )
  dimnames( lprof ) <- list( c( names( theta[1:3] ), "pp", "lprof" ), 1:nsupp )
  for( ii in 1:nsupp ) {
    theta[4] <- psupp[ii]

    ## evaluate profile loglikelihood
    thetas[1:4,ii] <- theta <- em.ggb( xx, yy, theta, theta[1:3], 0,
                                      printit = TRUE )
    lprof[5,ii] <- loglik( theta, xx, yy )
  }
  lprof
}
#######################################################################
s.marg <- function( xx, yy,
                   aa = 22.8, a0 = 1.08, nuA = .01, nu0 = .159, p = .064 )
{
  ## Compare empirical distribution of each color, say xx, or yy, against
  ## its fitted distribution

  ##     cy3/cy5     a    a0      nu.g   nu.r   nu     p
  ##MN1   1.27     32.9  1.33    0.011  0.016  0.233  0.033
  ##MN2a  1.27     22.8  1.08    0.010  0.014  0.159  0.064
  ##MN2b  1.30     15.1  0.84    0.009  0.008  0.174  0.050
  ##MN3a  1.64      3.9  1.90    9.15   4.12   1.29   0.212
  ##MN3b  1.60      2.5  1.93   18.2    6.38   2.36   0.343

  supp <- seq( min(x), max(x), length=500 )
  logmargA <- lgamma(aa+a0) - lgamma(aa) - lgamma(a0) +
    a0*log(nuA) + (aa-1)*log(supp) - (aa+a0)*log(supp+nuA)
  logmarg0 <- lgamma(aa+a0) - lgamma(aa) - lgamma(a0) +
    a0*log(nu0) + (aa-1)*log(supp) - (aa+a0)*log(supp+nu0)
  
  p * exp(logmargA) + (1-p) * exp(logmarg0)
}
#######################################################################
shrinkplot <- function( xx, yy, fits = s.fits( xx, yy ), chip="Control")
{
## Fits the Gamma-Gamma model (like s.five from earlier)

  xhat <- xx + fits[1]
  yhat <- yy + fits[1]
  plot( xx, yy, log="xy", pch=".", xlab="Cy3", ylab="Cy5",
       xlim=lims, ylim=lims )
  text( .01, 100, chip, cex=.8, adj=0 )
  
  mm <- length(xx)
  for( i in 1:mm )
    lines( c( xx[i], xhat[i] ), c(yy[i],yhat[i]), lwd=.2)
  invisible()
}
#######################################################################
oddsplot <- function( x, y, theta, by.level = 10,
                     rotate = FALSE, offset = 0,
                     main = "", xlab = xlabs, ylab = ylabs,
                     col = NULL, cex = c(.25,.75),
                     shrink = FALSE,
                     lims = range( c( x, y )))
{
  ## Plot odds curve for Gamma Gamma Bernoulli model

  ## truncate negative values for evaluation
  tmp <- x > -offset
  x <- x + offset
  x[!tmp] <- min( x[tmp] ) / 2
  tmpy <- y > -offset
  y <- y + offset
  y[!tmpy] <- min( y[tmpy] ) / 2
  tmp <- !( tmp & tmpy )
  if( any( tmp ))
    warning( paste( sum( tmp ), "probes truncated to", offset )) 
  rm( tmp )

  logbf <- lod.ggb(x,y,theta=theta)

  if( shrink ) {
    x <- x + theta[3]
    y <- y + theta[3]
  }
  ylabs <- "Cy5"
  xlabs <- "Cy3"  
  if( rotate ) {
    tmp <- sqrt( x * y )
    y <- y / x
    x <- tmp
    rm( tmp )
    ylabs <- paste( ylabs, xlabs, sep = " / " )
    xlabs <- "Average Intensity"
    if( missing( lims )) {
      xlim <- range( x )
      ylim <- range( y )
    }
  }
  else {
    xlim <- ylim <- lims
  }
  par( pty = "s" )

  plot( x[1], y[1], log="xy", xlab=xlab, ylab=ylab, xlim=xlim,
       ylim=ylim, type="n" )
  title( main )
  ##  usr <- par( "usr" )
  ##  text( 10^( usr[1]+ strwidth("abc") ), 10^((usr[3]+3*usr[4])/4), main,
  ##     cex=.8, adj=0 )
	
  ## report points with LOD > 0
  tmp <- logbf >= 0
    if( missing( col ) | is.null( col )) {
    col <- rep( "black", length( x ))
    col[tmp] <- "blue"
  }
  if( length( col ) != length( x ) & length( col ) != 1 )
    col <- col[1]
  for( i in unique( col )) {
    coli <- ( i == col ) & tmp
    if( any( coli ))
      points( x[coli], y[coli], cex=cex[2], col=i )
    coli <- ( i == col ) & !tmp
    if( any( coli ))
      points( x[coli], y[coli], cex=cex[1], col=i )
  }

  ## contour lines
  if( rotate ) {
    abline( h = 1, lty = 2, col = "red" )
    assign( "rlod.ggb", function( z, w, theta ) {
      w <- sqrt( w )
      lod.ggb( z / w, z * w, theta ) } )
    fun <- "rlod.ggb"
  }
  else {
    abline( 0, 1, lty = 2, col = "red" )
    fun <- "lod.ggb"
  }
  vec <- seq( log10( lims[1] ), log10( lims[2] ), length = 100 )

  bf <- if( shrink )
    outer( 10^vec - theta[3], 10^vec - theta[3], fun, theta = theta )
  else
    outer( 10^vec, 10^vec, fun, theta = theta)
  
  ## filled.contour(10^vec,10^vec,bf,levels=c(0,5),col=c("lightgray","white"),
  ##   save=TRUE, plotit=TRUE, add=TRUE, labex=0, lwd=2 )
  
  ## contours at 0,1,2 LOD
  contour(10^vec,10^vec,bf,levels=0,
          save=TRUE, plotit=TRUE, add=TRUE, labex=0, lwd=1, col = "red", lty = 3 )
  if( max( logbf ) >= 1 & by.level > 0 )
    contour(10^vec,10^vec,bf,
            levels=seq(0,floor(max(logbf)),by=log10(by.level))[-1],
            save=TRUE, plotit=TRUE, add=TRUE, labex=0, lwd=1, lty = 3 )
  
  ## box()
  ## tt <- x/y
  ## chat <- sqrt( mean( (tt-1)^2/(1+tt^2) ) )
  ## tmp01 <- chen.poly(chat,err=.01)
  ## abline( -log(tmp01[1]), 1, lty=2, lwd=1.5, err=(-1) )
  ## abline( -log(tmp01[2]), 1, lty=2, lwd=1.5, err=(-1) )
  
  invisible( logbf )
}
#######################################################################
lodprobes <- function( xx, yy, theta, lod, probes, col = 1, lowlod = 0,
                      offset = 0 )
{
  tmp <- xx > -offset & yy > -offset
  xx <- xx[tmp] + offset
  yy <- yy[tmp] + offset
  if( any( !tmp ))
    warning( paste( sum( !tmp ), "probes dropped with values below", offset ))
  rm( tmp )

  tmpc <- lod >= lowlod
  lod.order <- order( -lod[tmpc] )
  
  ## everything is ordered by LOD score

  ## probe names
  add.probes <- as.character( probes[tmpc] )
  add.probes <- add.probes[lod.order]
  
  ## LOD score
  lod.probes <- data.frame( probe = add.probes,
                           LOD = -sort(-lod)[seq(length(add.probes))] )

  ## ratio of xx to yy
  lod.probes$ratio <- c(( xx[tmpc] + theta[3] ) /
                        ( yy[tmpc] + theta[3] ))[lod.order]

  ## signed LOD score
  lod.probes$LOD[lod.probes$ratio<1] <- -lod.probes$LOD[lod.probes$ratio<1]

  ## inverse ratio
  lod.probes$inverse <- 1 / lod.probes$ratio

  ## round off numbers to 3 decimal places
  lod.probes[,-1] <- round(lod.probes[,-1],3)

  ## colors from plot
  if( length( col ) == length( xx ))
    lod.probes$col <- c(col[tmpc])[lod.order]

  lod.probes
}
#######################################################################
s.check0 <- function( xx, yy, theta1, theta2, chip = "Control" )
{
  lims <- c(.0065,1208)

  ## work it on the natural log scale
  supp <- seq( log(lims[1]), log(lims[2]), length=100 )

  lden <- function(x,aa,a0,nu) {
    ## returns log density of natural log of intensity
    lgamma(aa+a0) - lgamma(aa) - lgamma(a0) +
      a0*log(nu) + (aa-1)*log( exp(x) ) - 
        (aa+a0)*log( exp(x) +nu) + x 
  }

  hist( log( c(xx,yy) ), 50, prob=TRUE, ylim=c(0,.46), xlab="",ylab="",
       xlim=c(-5,8), cex=.9 )

  aa <- theta1[1]; a0 <- theta1[2]; nu <- theta1[3]
  logmarg <-  lden( supp, aa, a0, nu )
  lines( supp, exp(logmarg), lty=1 , lwd=2)

  aa <- theta2[1]; a0 <- theta2[2]; nu <- theta2[3]
  logmarg <-  lden( supp, aa, a0, nu )
  lines( supp, exp(logmarg), lty=2 , lwd=2)

  text( -5, .35, adj=0, chip, cex=.8 )
  invisible()
}
#######################################################################
s.check1 <- function( xx, yy, theta, chip = "Control" )
{
  ## Check the fit of the Gamma-Gamma-Bernoulli model by
  ## looking at (R-G)/(R+G) for spots deemed to not change.

  supp <- seq(.001,.999,length=150)

  logbf <- lod(xx,yy,theta=theta)
  ind <- (logbf < 0 )
  xx <- xx[ind]
  yy <- yy[ind]
  stat <- .5*( (xx-yy)/(xx+yy) + 1 )
  hist(stat,50,prob=TRUE,ylim=c(0,6.7), cex=.9 )

  den <- dbeta(supp,theta[1],theta[1])
  lines(supp,den,lwd=2)

  text(.1,5,chip,adj=0,cex=.8)
  invisible()
}
#######################################################################
s.check2 <- function( foo, xa, ya, thetaa, xb, yb, thetab,
                     spots = dimnames( foo)[[1]] )
{
  ## Look at the genes with high LOD compared to predicted changes
  ##from Craig's paper 

  fix.spots <- function( x, y, theta, spots ) {
    ## fix the edges
    x[x<0] <- 0
    y[y<0] <- 0
    bfa <- lod.ggb(x,y,theta=theta)

    ## skim the top
    inda <- (1:4290)[bfa>0]
    tmpa <- bfa[inda]
    orda <- order( -tmpa )
    ( spots[inda] )[orda]
  }
  spa <- fix.spots( xa, ya, thetaa, spots )
  spb <- fix.spots( xb, yb, thetab, spots )
  
  blah <- outer( spa, spb, "==" )
  bar <- apply(blah,1,any)  # the longer dimension
  spa[bar]
}

#######################################################################
## ftp://ftp.biostat.wisc.edu/pub/newton/npvolume/s.tack
## Beckett Diaconis Tack Data
## nsuccess <- c( rep(1,3), rep(2,13), rep(3,18), rep(4,48), 
##    rep(5,47), rep(6,67), rep(7,54), rep(8,51), rep(9,19) )
## ntrials <- 9; N <- length( nsuccess )
##   ## Binomial likelihood
## 	db2 <- function(y,prob,n){return(dbinom(y,n,prob))}
## 	lik <-  outer(nsuccess,grid,FUN="db2",n=ntrials) 
## gg = dbeta(grid,shape1=.5,shape2=.5),
gammaden <- function( x, a, b )
{
  b^a * x^(a-1) * exp( -x*b ) / gamma( a )
}
## recurbayes <- function( x, theta, domain = c(.01,.99), ngrid = 100
##   grid = seq(domain[1],domain[2],length=ngrid),
##   lik, gg = gammaden( exp(grid), theta[1], theta[3] ),
##   xlab="tack success probability", 
##   ylab="posterior predictive density" )
## {
## ## grid = support of mixing distribution 
## ## gg = prior guess
## ## lik = likelihood for data
##
##   N <- length( x )
##   alpha <- 1/3 
##
## ## A weight sequence
##   weight <- 1/sqrt((alpha+1)*(alpha+1:N))
##
## ## Process tacks in random order
##   ord <- sample( 1:N )
##   delta <- grid[2]-grid[1]
##
## ## Recursion yields approximate Bayes estimate gg
##   for( i in 1:N ) {
##     post <- lik[ord[i],]*gg
##     post <- ( post/sum(post) )/delta
##     gg <- gg*( 1-weight[i] ) + weight[i]*post 
##   }
## ## Good idea to repeat loop to see variation over orderings.
##
## ## Estimated predictive density
##   plot( grid, gg, type="l", xlab=xlab, ylab=ylab )
##
##   invisible( gg )
## }
#######################################################################
## This code uses a nonparametric Bayesian predictive recursion to
## estimate the mixing distribution of scale parameters for Gamma
## distributed array data.  See Newton and Zhang (1999) Biometrika, 86,
## 15-26 for more about this recursive algorithm, or go to
## www.stat.wisc.edu/~newton/
##
## The purpose of this calculation is to diagnose inadequacies of
## the Gamma mixing assumption in the Gamma-Gamma gene expression model.
########################################################################
predrecur <- function( xx, theta = c(32.9,1.33,0.01), gridlim = c(.0001,1) )
{
## This file takes in one set of array measurements in the vector xx
  N <- length( xx )

## theta: Observation component of the model is treated as known
## Gamma(aa,theta)
## Take as a prior guess of the mixing distribution for theta a Gamma(a0,nu)
## as estimated from the G-G model

## upp: Support of mixing distribution for random effects theta
## plug in here an upper support limit for the mixing distribution
  grid <- seq( gridlim[1], gridlim[2], length = 150 )
  delta <- grid[2] - grid[1]

## Gamma prior guess
  gg <- dgamma( grid, shape = theta[2], scale = ( 1 / theta[3] ))
  alpha <- 1
  g0 <- gg

## Gamma likelihood
  dg2 <- function( y, theta, shape )
    dgamma( y, shape = shape, scale = ( 1 / theta ))
  lik <-  outer( xx, grid, FUN = "dg2", shape = theta[1] )

## Recursion yields approximate Bayes estimate gg
  weight <- 1/sqrt( ( alpha + 1 ) * ( alpha + 1:N )) # A weight sequence
  ord <- sample( 1:N )    # Process genes in random order
  for( i in 1:N ) {
    post <- lik[ord[i],] * gg
    post <- ( post / sum(post) ) / delta
    gg <- gg * ( 1 - weight[i] ) + weight[i] * post 
  }
  plot( grid, gg, type="l", xlab="scale", 
       ylab="posterior predictive density" )
  lines( grid, g0, lty=2 )   # prior guess
  invisible( gg )
}
## (Repeat loop to see variation over orderings.)
## (We recommend averaging over a half dozen or so orderings---)

## Plot the estimated predictive density for the scale parameter theta
## of a future spot.


##########################################################################

