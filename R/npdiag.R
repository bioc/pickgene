###########################################################################
##
## s.npdiag: S code to do a nonparametric mixing diagnostic for microarrays
##
##     Copyright (C) 2000 Michael A. Newton.
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
#
# This Splus/R code uses a nonparametric Bayesian predictive recursion to
# estimate the mixing distribution of scale parameters for Gamma
# distributed array data.  
#
# The pupose of this calculation is to diagnose inadequacies of
# the Gamma mixing assumption in the Gamma-Gamma-Bernoulli
# gene expression model.  
#
# See Newton, Kendziorski, Richmond, Blattner, and Tsui 1999,
# http://www.stat.wisc.edu/~newton/papers/abstracts/btr139a.html
# where the hierarchical parametric models for gene expression data
# are presented.
#
# See Newton and Zhang (1999) Biometrika, 86, 15-26 for more about 
# the recursive algorithm used in the code below, or go to
# http://www.stat.wisc.edu/~newton/research/npb.html
#
#
########################################################################

npdiag <- function( xx, yy, aa, a0, nu, pp ) {
                                        # xx,yy raw expression measurement
                                        # vectors of equal length
                                        # Parameters of the fitted GGB model:
                                        # aa shape, observation component
                                        # a0 shape, inner component
                                        # nu scale, inner component
                                        # pp proportion of changed genes

  N <- length(xx)

                                        # grid range to approximate mixing
  ends <- log( (1/nu)*qgamma( c(.001,.999), shape=a0 ) )

                                        # Support of mixing distribution 
  grid <- seq(ends[1],ends[2],length=350)
  delta <- grid[2]-grid[1]

                                        # log-Gamma prior guess
  dgam <- function(y,shape,scale){ return( scale*dgamma(y*scale,shape) ) }
  g0 <- exp(grid)*dgam(exp(grid),shape=a0,scale=nu )
  gg <- g0

                                        # Gamma likelihood with logged scale
  dg2 <- function( y,theta,shape ) {
    return( dgam( y, shape = shape, scale = exp( theta )))
  }

                                        # Recursion yields approximate
                                        # Bayes estimate gg
  alpha <- 1
  weight <- 1/sqrt((alpha+1)*(alpha+1:N)) # A weight sequence

  ord <- sample( 1:N )    # Process genes in random order
  for( i in 1:N ) {
    j <- ord[i]
    z <- c( xx[j], yy[j] )
                                        # Joint prob( data_j | mixing parameter )
    lik <-  z*outer(z,grid,FUN="dg2",shape=aa) 
    tmp <- c( lik[1,] * lik[2,] )
    post <- tmp*gg
    post <- ( post/sum(post) )/delta
    gg <- gg*( 1-weight[i] ) + weight[i]*post 
    print(i)
  }
                                        # Repeat loop to see variation over orderings.


                                        # Take a look at the results

  plot( grid, gg, type="l" )   # nonparametric estimate
  lines( grid, g0, lty=2 )     # parametric estimate/prior guess
}

