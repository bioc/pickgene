#####################################################################
##
## $Id: yandell.R,v 1.2 2003/12/30 19:17:27 jgentry Exp $
##
##     Copyright (C) 2001 Brian S. Yandell
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
rangene <- function (n = 10000, center = 4, spread = 2, contamination = 0.05,
    alpha = c(1,1), noise = 0.5, omega = c(20,20))
{
    contam <- round(n * contamination)
    data <- list()
    r <- rnorm(n, center, spread)
    rr <- rank(r) / ( 1 + n )
    ## contaminate preferentially genes with smaller intensities
    cc <- seq( n )[ contam >= rank( runif( n )^2 * (2 - rr )) ]
    ## half of contaminated genes upregulate
    upreg <- .5 < runif( contam )
    for (i in 1:2) {
        data[[i]] <- r
        s <- cc[upreg]
        if( length( s ))
          data[[i]][s] <- data[[i]][s] + rnorm(length(s),3-rr[s],0.25)
        upreg <- !upreg
    }
    for (i in 1:2)
        data[[i]] <- data[[i]] + rnorm(n, 0, noise)

    known <- data <- as.data.frame(data)
    for (i in 1:2) {
        data[[i]] <- alpha[i] * exp(data[[i]])
        data[[i]] <- data[[i]] + rnorm(n, 0, omega[i])
    }
    list(observed = data, true = known,contam=cc)
}
### Find top ranked genes
toprankgene <- function( yy, n = 500 )
{
  tmpy2 <- rank( apply( yy$true, 1, function(x) -abs( diff( x )) ))
  tmpy <- pickgene(yy$obs, npickgene = n )$pick[[1]]$probe

  count <- rep(NA,n)
  for( i in 1:n ) {
    count[i] <- sum( tmpy2[ tmpy[1:i] ] <= i )
  }
  count
}
### Yi Lin's original random gene generator
### Classical contamination model for robust statistics
orangene <- function( n = 10000, center = 4, spread = 2, contamination = .05,
         alpha = c(3,2), noise = 0.1, omega = c(20,30) )
{
  contam <- round( n * contamination )

  data <- list()
  r <- rnorm(n,center,spread)
  for( i in 1:2 ) {
    ## gene expression
    data[[i]] <- r

    ## contamination = differential expression
    data[[i]][seq(contam)] <- data[[i]][seq(contam)] + rnorm(contam)

    ## intrinsic noise
    data[[i]] <- data[[i]] + rnorm(n,0,noise)
  }
  known <- data <- as.data.frame( data )
  for( i in 1:2 ) {
    ## attenuation
    data[[i]] <- alpha[i] * exp( data[[i]] )

    ## measurement error
    data[[i]] <- data[[i]] + rnorm(n,0,omega[i])
  }
  list( observed = data, true = known )
}
### Organize scores from pickgene structure
pickedhist <- function( pick, show = names( pick ), title = show, p1 = .05,
                       plotit = TRUE, rotate = FALSE, mfrow = c(nr,nc), bw = NULL )
{
  pick <- pick$score
  if( is.numeric( show ))
    show <- names( pick )[show]
  else {
    if( any( is.na( match( show, names( pick )))))
      stop( paste( "show choices do not match contrasts:",
                  paste( names( pick ), collapse = ", " )))
  }
  n <- nrow( pick[[1]] )
  if( plotit ) {
    nr <- min( 2, ceiling( length( show ) / 3 ))
    nc <- min( 3, ceiling( length( show ) / nr ))
    if( !is.null( mfrow ))
      par( mfrow = mfrow, pty = "s" )
    tmplvl <- qnorm( adjustlevel( n * 2, .05 ) / 2, lower.tail
                    = FALSE)
    names( title ) <- show
    if( !is.null( bw )) {
      bw <- array( bw, length( show ))
      names( bw ) <- show
    }
    for( i in show ) {
      if( is.null( bw ))
        tmp <- density( pick[[i]]$score )
      else
        tmp <- density( pick[[i]]$score, bw = bw[i] )
      cat( "density bandwidth:", tmp$bw, "\n" )
      if( rotate )
        plot( tmp$y, tmp$x, type = "l", ylab = title[i],
             xlab = "Relative Frequency", main = "" )
      else
        plot( tmp$x, tmp$y, type = "l", xlab = title[i],
             ylab = "Relative Frequency", main = "" )
      tmpd <- dnorm( tmp$x )
      if( rotate )
        lines( tmpd, tmp$x, lty = 2, col = "blue" )
      else
        lines( tmp$x, tmpd, lty = 2, col = "blue" )
      tmphi <- 2 * tmpd < tmp$y
      xx <-  range( tmp$x[ abs( tmp$x ) < tmplvl & !tmphi] )
      if( rotate ) {
        lines( tmp$y * tmphi, tmp$x, lty = 4, col = "red" )
        abline( h = c(-1,1) * tmplvl, lty = 3, col = "red" )
        axis( 2, xx, labels = round( xx, 1 ))
        lines(( tmp$y - (1 - p1) * tmpd ) / p1, tmp$x, lty = 5, col = "purple" )
      }
      else {
        lines( tmp$x, tmp$y * tmphi, lty = 4, col = "red" )
        abline( v = c(-1,1) * tmplvl, lty = 3, col = "red" )
        axis( 1, xx, labels = round( xx, 1 ))
        lines( tmp$x, ( tmp$y - (1 - p1) * tmpd ) / p1, lty = 5, col = "purple" )
      }
    }
  }
  invisible( tmp )
}
### Organize scores from pickgene structure
pickedchisq <- function( pick, show = names( pick ),
                        title = "Squared Distance",
                        plotit = TRUE, alpha = .05 )
{
  pick <- pick$score
  if( is.numeric( show ))
    show <- names( pick )[show]
  else {
    if( any( is.na( match( show, names( pick )))))
      stop( paste( "show choices do not match contrasts:",
                  paste( names( pick ), collapse = ", " )))
  }
  n <- nrow( pick[[1]] )
  if( plotit ) {
    tmplvl <- sqrt( qchisq(( 1 - alpha ) ^ ( 1 /n ), length( show )))
    sumsq <- 0
    for( i in show )
      sumsq <- sumsq + pick[[i]]$score^2
    tmp <- density( sumsq )
    plot( tmp$x, tmp$y, type = "l", log = "x", xlab = title,
             ylab = "Relative Frequency", main = "" )
    tmpd <- dchisq( tmp$x, length( show ))
    lines( tmp$x, tmpd, col = "blue" )
    tmphi <- 2 * tmpd < tmp$y
    xx <-  max( tmp$x[ tmp$x < tmplvl & !tmphi ] )
    lines( tmp$x, tmp$y * tmphi, col = "red" )
    abline( v = tmplvl, col = "red" )
    axis( 1, xx, labels = round( xx, 1 ))
    lines( tmp$x, tmp$y - tmpd, col = "purple" )
  }
  invisible( tmp )
}
holms <- function( x, alpha = .05, cut = TRUE ) {
  n <- length( x )
  x <- sort( 2 * ( 1 - pnorm( abs( x ))))
  x[x==0] <- min( x[x>0] ) / 2
  if( cut )
    x <- x[ x <= alpha ]
  nx <- length( x )
  plot( x, log = "y", xlab = "rank order", ylab = "raw p-value" )
  axis( 2, alpha )
  abline( h = c(alpha,alpha/n), col = "red" )
  text( nx, 0.5 * ( alpha / n ), "Bonferroni", adj = 1, col = "red" )
  abline( h = 1-(1-alpha)^(1/n), col = "blue" )
  text( nx, 2 * ( 1-(1-alpha)^(1/n)), "Sidak", adj = 1, col = "blue" )
  lines( alpha / ( 1 + n - seq( x )), col = "purple" )
  text( nx, alpha / sqrt(n), "Holms", adj = 1, col = "purple" )
  abline( v = max( seq( x )[ x < alpha / n ] ), col = "blue", lty = 2  )
  abline( v = max( seq( x )[ x < alpha / ( 1 + n - seq( x )) ] ),
    col = "purple", lty = 2 )
  points( seq( x ) / ( 1 + n ), col = "grey" )
  invisible( x )
}
### Organize scores from pickgene structure
pickedscore <- function( pick, description, show = 1:2, alpha = .05,
                        xlab = show[1], ylab = show[2], main = "",
                        mfrow = c(1,1))
{
  score <- list()
  for( i in names( pick$score ))
    score[[i]] <- pick$score[[i]]$score
  score <- data.frame( score )

  show <- array( show, 2 )
  if( is.numeric( show ))
    show <- names( pick$pick )[show]
  else {
    if( any( is.na( match( show, names( pick$pick )))))
      stop( paste( "show choices do not match contrasts:",
                  paste( names( pick$pick ), collapse = ", " )))
  }
  if( !is.null( mfrow ))
    par( mfrow = mfrow, pty = "s" )

  library(MASS)
  eqscplot( score[[show[1]]], score[[show[2]]], type = "n",
           xlab = xlab, ylab = ylab )
  title( main )
  points( score[[show[1]]], score[[show[2]]], cex = .25 )
  n <- length( score[[1]] )
  tmplvl <- qnorm( adjustlevel( n * 2, alpha ) / 2, lower.tail
                  = FALSE)
  abline( h = c(-1,1)*tmplvl, v = c(-1,1) * tmplvl, col = "red" )
  tmp <- sqrt( qchisq(( 1 - alpha ) ^ ( 1 /n ), 2 ))
  x <- seq( -tmp, tmp, length = 200)
  y <- sqrt( tmp^2 - x^2 )
  lines(x,y,col="blue")
  lines(x,-y,col="blue")
  tmp <- abs( score[[show[1]]] ) > tmplvl | abs( score[[show[2]]] ) >
    tmplvl
  points( score[[show[1]]][tmp], score[[show[2]]][tmp], cex = .75, col =
         "blue" )

  probes <- character()
  for( i in names( pick$pick ))
    probes <- c( probes, as.character( pick$pick[[i]]$probe ))
  probes <- sort( unique( probes ))
  score <- score[ match( probes, pick$score[[1]]$probe, nomatch = 0 ), ]
  dimnames( score ) <- list( probes, names( score ))

  pick <- pick$pick
  fold <- pvalue <- matrix( NA, length( probes ), length( pick ),
                           dimnames = list( probes, names( pick )))
  for( i in names( pick )) {
    fold[ as.character( pick[[i]]$probe ), i ] <- pick[[i]]$fold
    pvalue[ as.character( pick[[i]]$probe ), i ] <- pick[[i]]$pvalue
  }
  ll <- list( score = score, fold = fold, pvalue = pvalue )
  ## add descriptions if provided
  if( !missing( description ))
    ll$description <- description[probes]
  ### want to add probe description to list
  ll
}
pickedpair <- function( x, columns, description,
                       probe = "Probe.Set",
                       renorm = c(sqrt(2),sqrt(6)),
                       pick = pickgene( x[,columns],
                         x[,probe], ...,
                         renorm = renorm, plotit = FALSE ),
                       main = "", ... )
{
  par(mfrow=c(2,2),pty = "s" )
  score <- pickedscore( pick, description, mfrow = NULL,
                              xlab = "Additive Effect",
                              ylab = "Dominance Effect",
                              main = main )
  pickedhist( pick, 2, title = "Dominance Effect", mfrow = NULL,
             rotate = TRUE )
  pickedhist( pick, 1, title = "Additive Effect", mfrow = NULL )
  invisible( score )
}
robustbox<-function(y,x,nslice=400,
  xlab = "Log Average Intensity", ylab = "Standardized Difference",
  shrink = FALSE,
  crit = qnorm( adjustlevel( n , overalllevel ) / 2, lower.tail = FALSE),
  overalllevel = .05, cex = 0.1, lwd = 2, plotit = TRUE )
{
  n<-length(x)
  ylim <- if( shrink )
    c(-1,1) * crit
  else
    range( y )

  ox<-order(x)
  x<-x[ox]
  y<-y[ox]
  rx<-rank(x)
  k<-floor(n/nslice)
  slicef<-factor(cut(rx, breaks = c(k*(0:(nslice-1)),n)))
  slicex<-unlist(tapply(x, slicef, median), recursive = TRUE,  use.names=FALSE)

  slicebox <- t( matrix( unlist( tapply( y, slicef, function(x)
    boxplot.stats(x)$stats ), recursive = TRUE, use.names=FALSE), 5 ))
  robmed <- apply( slicebox, 2, median )
  plot(x,y, type = "n", xlab = xlab, ylab = ylab, ylim = ylim )
  if( plotit )
    points(x,y,cex = cex)
  for( i in seq( ncol( slicebox ))) {
    if( diff( range( slicebox[,i] )) > 0 ) {
      tmp <- smooth.spline( slicex, slicebox[,i] )
      tmp <- predict( tmp, x )$y
      lines( x, tmp, col = "blue", lwd = lwd )
    }
    else
      lines( range( x ), rep( robmed[i], 2 ), col = "red", lwd = lwd )
  }
  invisible()
}
