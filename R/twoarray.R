#####################################################################
##
## $Id$
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
##
# Using the normalization from Richmond et al. (1999)
# File to read data for subsequent analysis
# This one normalizes to average intensity first

twoarray.norm <- function( foo, ..., conditions = c("Cy3","Cy5"),
                       reduce = FALSE, identifier = "identifier" ) {

  foonames <- list( ... )
  if( !is.null( names( foonames )) & missing( conditions ))
    conditions <- names( foonames )
  xnorm <- foonames[[1]]
  ynorm <- foonames[[2]]

  ## Background adjustment (very simple)
  x <- foo[[xnorm[1]]] - foo[[xnorm[2]]]
  y <- foo[[ynorm[1]]] - foo[[ynorm[2]]]

  ## Normalization
  ## Rescale to help with underflow problem 10^5 (does not affect shape params)
  x <- 100000 * x / sum( x[x>0] )
  y <- 100000 * y / sum( y[y>0] )

  if( reduce ) {
    ok <- x > 0 & y > 0
    l <- data.frame( foo[[1]][ok], x[ok], y[ok] )
  }
  else {
    l <- data.frame( foo[[1]], x, y )
  }
  names( l ) <- c( identifier, conditions )
  l
}

## compare oddsplot and pickgene
twoarray.plot <- function( mydata,
                          main = deparse( substitute( mydata )),
                          theta,
                          conditions = c("Cy3","Cy5"),
                          identifier = "identifier" ) {

  par(mfrow=c(2,2))

  ## oddsplot: Newton et al. (2001)
  hist( log( mydata[[conditions[1]]] )[mydata[[conditions[1]]]>0],
       breaks = 100, xlab = "Log Adjusted Condition 1")
  hist( log( mydata[[conditions[2]]] )[mydata[[conditions[2]]]>0],
       breaks = 100, xlab = "Log Adjusted Condition 2")
  plot( log( mydata[[conditions[1]]] ), log( mydata[[conditions[2]]] ),
       xlab = paste( "Log Adjusted", conditions[1] ),
       ylab = paste( "Log Adjusted", conditions[1] ), cex = 0.4 )
  title( main )
  xxx <- qnorm( rank( mydata[[conditions[1]]] ) /
               ( 1 + length( mydata[[conditions[1]]] )))
  yyy <- qnorm( rank( mydata[[conditions[2]]] ) /
               ( 1 + length( mydata[[conditions[2]]] )))

  if( !missing( theta ))
    do.oddsplot( mydata, main, rotate = TRUE, theta = theta )
  else
    do.oddsplot( mydata, main, rotate = TRUE )

  par( mfcol = c(3,2))

  ## pickgene: Lin et al. (2001)
  plot( log( mydata[[conditions[1]]] ), log( mydata[[conditions[2]]] ),
       xlab = paste( "Log Adjusted", conditions[1] ),
       ylab = paste( "Log Adjusted", conditions[2] ), cex = 0.4 )
  title( paste( "Log of", main ))
  abline( 0, 1, col = "red", lty = 2 )
  tmp <- mydata[[conditions[1]]] > 0 & mydata[[conditions[2]]] > 0
  a <- robustscale(( log( mydata[[conditions[1]]] )
                    - log( mydata[[conditions[2]]] ))[tmp],
                   ( log( mydata[[conditions[1]]] )
                    + log( mydata[[conditions[2]]] ))[tmp] )
  qqnorm(( a$y - a$center ) / a$scale, cex = 0.4 )
  abline( 0, 1, col = "red", lty = 2 )

  pickgene( mydata[,conditions], mydata[[identifier]], mfrow = NULL )
  plot(xxx, yyy,
       xlab = paste( "Rank Adjusted", conditions[1] ),
       ylab = paste( "Rank Adjusted", conditions[2] ) ,cex=0.4)
  title( paste( "Rank of", main ))
  abline( 0, 1, col = "red", lty = 2 )
  a <- robustscale( xxx-yyy, xxx+yyy )
  ## hist(( a$y - a$center ) / a$scale, breaks = 100 )
  qqnorm(( a$y - a$center ) / a$scale, cex = 0.4 )
  abline( 0, 1, col = "red", lty = 2 )
  pickgene( mydata[,conditions], mydata[[identifier]],
    rankbased = TRUE, mfrow = NULL )

}
do.oddsplot <- function(data,
                        main = substitute( data ),
                        theta = c(2,2,2,.4),
                        col = NULL,
                        redo = missing( theta ),
                        conditions = c("Cy3","Cy5"),
                        identifier = "identifier", ... ) {
  if( redo )
    theta <- em.ggb(data[[conditions[1]]], data[[conditions[2]]],
                    theta, theta[1:3], print = TRUE )
  lod <- oddsplot( data[[conditions[1]]], data[[conditions[2]]], theta,
                  xlab = conditions[1], ylab = conditions[2],
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
