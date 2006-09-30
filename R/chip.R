#####################################################################
##
## $Id$
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
makecont <- function( x, y, size = 41, cex = .1,
                     levels = c(1,5,10,50))
{
  isna <- is.na(x) | is.na(y)
  if( all( isna ))
    stop("all data missing")
  x <- x[!isna]
  y <- y[!isna]
  rx <- range( x )
  ry <- range( y )
  z <- matrix( 0, size, size )
  xx <- 1 + floor( ( size + 1 ) * ( x - rx[1] ) / diff( rx ))
  xx[xx>size] <- size
  yy <- 1 + floor( ( size + 1 ) * ( y - ry[1] ) / diff( ry ))
  yy[yy>size] <- size
  for( i in seq( length( x ))) {
    z[xx[i],yy[i]] <- 1 + z[xx[i],yy[i]]
  }
  z <- 100 * z / length( x )
  oz <- order( z )
  z[oz] <- cumsum( z[oz] )
  oz <- z[ xx + size * ( yy - 1 ) ] < 1
  
  sx <- seq( rx[1], rx[2], length = size )
  sy <- seq( ry[1], ry[2], length = size )
  contour( sx,sy,z, levels = levels )
  
  points( x[oz], y[oz], cex = cex )
  invisible( z )
}
lod.plot <- function( data, x, y, theta,
                     filename = deparse( substitute( theta )),
                     probe = "Probe.Set",
                     xlab = x, ylab = y, ps = TRUE,
                     col = rep( "black", length( x )),
                     lowlod = 0,
                     ... )
{
  col <- as.character( col )
  if( ps )
    postscript( paste( filename, "ps", sep = "." ), horizontal = FALSE )
  par( omi = rep(.5,4))
  sink( paste( filename, "lod", sep = "." ))
  print( lodprobes( exp( data[[x]] ), exp( data[[y]] ), theta,
                   oddsplot( exp( data[[x]] ), exp( data[[y]] ),
                            theta, col = col, xlab=xlab, ylab=ylab,
                            chip=""),
                   data[[probe]], col, lowlod = lowlod ))
  sink()
  if( ps )
    graphics.off()
  invisible()
}
chipnorm <- function( xx, chip = rep( 1, length( xx )))
{
  chip <- as.factor( chip )
  chipmean <- tapply( xx, chip, mean, na.rm = TRUE )
  print( chipmean )
  for( i in levels( chip )) {
    tmp <- chip == i
    xx[tmp] <- xx[tmp] - chipmean[i]
  }
  list( xx = xx, mean = chipmean )
}
