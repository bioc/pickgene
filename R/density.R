#####################################################################
##
## $Id: density.R,v 1.1 2003/11/25 14:50:35 jgentry Exp $
##
##     Copyright (C) 2000 Brian S. Yandell
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
sixden <- function( x, y, align=FALSE, crit = 5, xlim=range(x-y, na.rm =TRUE),
                   dolog = TRUE, dif = x - y, ave = (x + y) / 2 )
{
  if( dolog ) {
    x <- log( x )
    y <- log( y )
    x[is.na(x)] <- 0
    y[is.na(y)] <- 0
  }
  par(mfrow=c(3,2))
  lx6 <- ceiling( length( dif ) / 6 )
  ranks <- order( order( ave ))
  for(i in seq( 0, length(dif), by = lx6)) {
    tmp <- ranks > i & ranks < i + lx6 + 1
    tmpd <- dif[tmp]
    tmpd <- tmpd[!is.na(tmpd)]
    ##		hist(tmpd, nc=200,xlab="log diff",freq=FALSE,xlim=xlim,
    tmpdd <- density( tmpd )
		
    if( align )
      plot( tmpdd$x, tmpdd$y, type="l",xlab="log diff", xlim=xlim,
           ylab="relative frequency")
    else
      plot( tmpdd$x, tmpdd$y, type="l",xlab="log diff", #xlim=xlim,
           ylab="relative frequency")
    title( main=paste( sum(tmp), "genes from",
             round(ave[ranks==i+1],2),
             "to", round(ave[ranks==i+lx6],2) ))
    tmpd <- tmpd[abs(tmpd-mean(tmpd,na.rm=TRUE))<5*sqrt(var(tmpd,na.rm=TRUE))]
    tmps <- c(mean(tmpd,na.rm=TRUE), sqrt(var(tmpd,na.rm=TRUE)))
    print(c(i,tmps))
    lines( sort(tmpd), dnorm(sort(tmpd),tmps[1],tmps[2] ),
          col="blue", lty = 2 )
    lines( rep( tmps[1],2),
          c(0,dnorm(tmps[1],tmps[1], tmps[2])),
          col = "red", lty = 2  )
    lines( rep( tmps[1]+tmps[2],2),
          c(0,dnorm(tmps[1]+tmps[2],tmps[1], tmps[2])),
          col = "red", lty = 2  )
    lines( rep( tmps[1]-tmps[2],2),
          c(0,dnorm(tmps[1]-tmps[2],tmps[1], tmps[2])),
          col = "red", lty = 2  )
  }
  par( mfrow = c(1,1))
}
denlines <- function( x, y, align=FALSE, crit = 5, xlim=range(x-y, na.rm =TRUE),
                     ylim=c(0,2.5), dolog = TRUE, dif = x - y, ave = (x + y) / 2,
                     numlines = 6, offset = 0 )
{
  if( dolog ) {
    x <- log( x )
    y <- log( y )
    x[is.na(x)] <- 0
    y[is.na(y)] <- 0
  }
  lx6 <- ceiling( length( dif ) / numlines )
  ranks <- order( order( ave ))
  plot( xlim, ylim, type="n",xlab="log diff",
       ylab="relative frequency")
  bump <- 0
  for(i in seq( 0, length(dif), by = lx6))
    {
      tmp <- ranks > i & ranks < i + lx6 + 1
      tmpd <- dif[tmp]
      tmpd <- tmpd[!is.na(tmpd)]
      tmpdd <- density( tmpd )
      bump <- bump + offset
      lines( tmpdd$x, tmpdd$y + bump)
    }
}
dencont <- function( x, y, align=FALSE, crit = 5, xlim=range(x-y, na.rm =TRUE),
                    ylim=c(0,2.5), dolog = TRUE, byranks = TRUE,
                    dif = x - y, ave = (x + y) / 2,
                    numlines = round( length( ave ) / 200 ),
                    levels.z = pretty( range( z ), 10 ))
{
  if( dolog ) {
    x <- log10( x )
    y <- log10( y )
    x[is.na(x)] <- 0
    y[is.na(y)] <- 0
  }
  lx6 <- ceiling( length( dif ) / numlines )
  dif <- dif[ order( ave )]
  ave <- sort( ave )
  rdif <- range( dif, na.rm = TRUE )
  ranks <- seq( length( ave ))

  z <- matrix( 0, 512, numlines )
  splits <- seq( 0, length(dif), by = lx6)
  j <- 0
  for(i in splits ) {
    tmp <- ranks > i & ranks < i + lx6 + 1
    tmpd <- dif[tmp]
    tmpd <- tmpd[!is.na(tmpd)]
    tmpdd <- density( tmpd, from = rdif[1], to = rdif[2] )
    j <- j + 1
    z[,j] <- 100 * tmpdd$y
  }
  if( byranks ) {
    tmpar <- par( yaxt = "n" )
    ytmp <- seq( 1, length( ave ), len = numlines )
  }
  else
    ytmp <- seq( min( ave, na.rm = TRUE ), max( ave, na.rm = TRUE ),
                len = numlines )
  xtmp <- seq( rdif[1], rdif[2], len = 512 )
  contour( xtmp, ytmp, z, levels = levels.z )
  if( byranks ) {
    par( yaxt = "s" )
    axis( 2, at = 1 + splits, labels = FALSE )
    tmp <- c( 1 + splits[ seq( 1, length( splits ),
                           by = floor( length( splits ) / 5 )) ],
             length( ave ))
    axis( 2, at = tmp, labels = round( ave[tmp], 2 ))
  }
  lines( apply( z, 2, function(z,x)x[z==max(z)][1], xtmp ),
        ytmp, lty = 2, col = "blue" ) 
  invisible( levels.z )
}
dencum <- function( x, y, align=FALSE, crit = 5, xlim=xlims,
                   ylim = ylims, dolog = TRUE, byranks = TRUE, standardize = FALSE,
                   dif = x - y, ave = (x + y) / 2, splineit = FALSE,
                   numlines = round( length( ave ) / 200 ), show = xlim,
                   levels.z = c(1,5,10,25,50,75,90,95,99))
{
  if( dolog ) {
    x <- log10( x )
    y <- log10( y )
    x[is.na(x)] <- 0
    y[is.na(y)] <- 0
  }
  probe <- seq( length( dif ))
  isna <- is.na( dif ) | is.na( ave )
  dif <- dif[!isna]
  ave <- ave[!isna]
  probe <- probe[!isna]

  lx6 <- ceiling( length( dif ) / numlines )
  dif <- dif[ order( ave )]
  probe <- probe[ order( ave )]
  ave <- sort( ave )
  rdif <- range( dif[ dif < Inf & dif > -Inf ], na.rm = TRUE )
  if( standardize )
    xlims <- c(-5,5)
  else
    xlims <- rdif
  ranks <- seq( length( ave ))

  cat( "Data divided into", numlines, "groups to determine quantiles\n" )
               
  z <- matrix( 0, length(levels.z), numlines )
  dimnames( z ) <- list( as.character( levels.z ), NULL )
  if( byranks ) {
    tmpar <- par( yaxt = "n" )
    ylims <- range( ranks )
  }
  else
    ylims <- range(ave, na.rm = TRUE)
  plot( xlim, ylim, type="n",xlab="log diff",
       ylab="log ave")

  splits <- seq( 0, length(dif), by = lx6)
  if( byranks ) {
    par( yaxt = "s" )
    axis( 2, at = 1 + splits, labels = FALSE )
    tmp <- c( 1 + splits[ seq( 1, length( splits ),
                           by = floor( length( splits ) / 5 )) ],
             length( ave ))
    axis( 2, at = tmp, labels = round( ave[tmp], 2 ))
  }
  nz <- nrow( z )
  if( byranks )
    ytmp <- seq( 1, length( ave ), len = 2 * numlines + 1 )
  else
    ytmp <- seq( min( ave, na.rm = TRUE ), max( ave, na.rm = TRUE ),
                len = 2 * numlines + 1 )

  showvals <- NULL
  j <- 0
  for( i in splits ) {
    tmp <- ranks > i & ranks < i + lx6 + 1
    tmp[tmp] <- !is.na( dif[tmp] )
    tmpd <- sort( dif[tmp] )
    j <- j + 1
    if( !byranks ) {
      ytmp[2*j] <- mean( ave[tmp] )
      ytmp[2*j+1] <- max( ave[tmp] )
    }
    z[,j] <- tmpd[ round( 0.01 * levels.z * length( tmpd )) ]
    if( standardize ) {
      mtmp <- 0; vtmp <- 1
      tmp2 <- tmp
      tmp2[tmp] <- dif[tmp] >= z[1,j] & dif[tmp] <= z[nz,j]
      if( any( tmp2 )) {
        mtmp <- mean( dif[tmp2] )
        vtmp <- var( dif[tmp2] )
        if( vtmp == 0 | is.na( vtmp ))
          vtmp <- 1
        else
          vtmp <- sqrt( vtmp )
      }
    }
    tmp[tmp] <- dif[tmp] < z[1,j] | dif[tmp] > z[nz,j]
    diftmp <- dif[tmp]
    if( standardize ) {
      z[,j] <- ( z[,j] - mtmp ) / vtmp
      if( length( diftmp )) {
        diftmp <- ( diftmp - mtmp ) / vtmp
        mtmp <-  diftmp > show[2] | diftmp < show[1]
        if( any( mtmp )) {
          showvals <- rbind( showvals, cbind( group = rep( j, sum( mtmp )),
            probe = ( probe[tmp] )[mtmp],
            ave = ( ave[tmp] )[mtmp],
            diff = diftmp[mtmp] ))
          diftmp[mtmp] <- pmin( xlim[2], pmax( xlim[1], diftmp[mtmp] ))
        }
      }
    }
    if( byranks )
      points( diftmp, ranks[tmp], cex = 0.5, col = "blue" )
    else
      points( diftmp, ave[tmp], cex = 0.5, col = "blue" )
  }
  if( length( showvals )) {
    cat( "Extreme values by group:\n")
    row.names( showvals ) <- seq( nrow( showvals ))
    print( showvals )
  }
  for( i in as.character( levels.z )) {
    tmpd <- c(rbind(z[i,],z[i,]))
    tmpd <- c(tmpd[1],(tmpd[-1]+tmpd[-2*numlines])/2,
              tmpd[2*numlines])
    lines( tmpd, ytmp )
    if( splineit ) {
      tmp <- spline( ytmp, tmpd )
      lines( tmp$y, tmp$x )
    }
    text( z[i,c(1,numlines)], ytmp[c(1,1+2*numlines)], i, cex = 0.5 )
  }
  invisible( showvals )
}
