###########################################################################
##
## $Id$
##
##     Copyright (C) 2000 Yi Lin and Brian S. Yandell
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
############################################################################
pickgene2 <- function( ... )
  pickgene.two( ... )
multipickgene<-function( ... )
  pickgene.poly( ... )
twowayanovapickgene <- function( x, fac1level, fac2level, ... )
  pickgene( x, faclevel = c(fac1level,fac2level), ... )
###########################################################################
## the following function robustly estimate the center and scale.
robustscale<-function(y,x,nslice=400,corcenter=TRUE,decrease=TRUE){
  n<-length(x)
  ox<-order(x)
  x<-x[ox]
  y<-y[ox]
  rx<-rank(x)
  k<-floor(n/nslice)
  slicef<-factor(cut(rx, breaks = c(k*(0:(nslice-1)),n)))
  slicex<-unlist(tapply(x, slicef, median), recursive = TRUE, use.names=FALSE)

  if (corcenter) {
    slicemedian<-unlist(tapply(y, slicef, median), recursive = TRUE, use.names=FALSE)
    slicemad<-unlist(tapply(y, slicef, mad), recursive = TRUE, use.names=FALSE)
  }
  else {
    slicemad<-unlist(tapply(y, slicef, mad, center=0), recursive = TRUE,
                     use.names=FALSE)
  }
  tmp <- slicemad > 0
  if( sum( tmp ) > 3 ) {
    preroblogscale <- smooth.spline(slicex[tmp],log(slicemad[tmp]))
    robscale <- exp(predict(preroblogscale,x)$y)
  }
  else {
    preroblogscale <- list( y = rep( log( median( slicemad[tmp] )),
                              length( slicemad )))
    robscale <- rep( median( slicemad ), length( x ))
  }
  if (decrease) {
    robscale[(ceiling(n/2)):1] <- cummax(robscale[(ceiling(n/2)):1])
    robscale[(ceiling(n/2)):n] <- cummin(robscale[(ceiling(n/2)):n])
  }
  if (corcenter) {
    if( length( slicex ) > 3 ) {
      prerobcenter <- smooth.spline(slicex,slicemedian/(exp(preroblogscale$y)))
      robcenter <- predict(prerobcenter,x)$y * robscale
    }
    else
      robcenter <- rep( median( slicemedian ), length( x ))
    adjustedrobscale <- mad((y-robcenter)/(robscale)) * robscale
  }
  else {
    robcenter <- rep(0,n)
    adjustedrobscale<-mad((y/robscale), center = 0) * robscale
  }
  return(list(center=robcenter,scale=adjustedrobscale,x=x,y=y))
}

###########################################################################
## significance level adjustment for multiple tests. (similar to Bonferroni)
adjustlevel<-function(ntest,alpha){
  singlelevel<- 1 - exp((log(1-alpha))/ntest)
  return(singlelevel)
}

###########################################################################
## the following function picks genes in 2 condition case with logratio (y)
## and logofproduct (intensity). It also serves as the building block for
## the later functions for multi-condition and ANOVA.
pickgene.two <- function ( y, intensity,
                          geneid = 1:length(y),
                          ## set elsewhere
                          singlelevel=0.0001,
                          ## automatic pick genes (or pick npickgene if >0)
                          npickgene = -1, meanrank = FALSE,
                          xlab="Average Intensity",
                          ylab="Trend", main = "",
                          plotit = TRUE, col = "blue",
                          negative = numeric( 0 ),
                          ... ) {
  n <- length( y ) - length( negative )
  if( length( negative )) {
    sc <- robustscale( y[-negative], intensity[-negative], ... )
    geneid <- geneid[ order( intensity[-negative] ) ]
  }
  else {
    sc <- robustscale( y, intensity, ... )
    geneid <- geneid[ order( intensity ) ]
  }
  z <- ( sc$y - sc$center ) / ( sc$scale )
  if( meanrank )
    logxy <- "y"
  else {
    logxy <- "xy"
    sc$x <- exp( sc$x )
  }
  sc$y <- exp( sc$y )
  zcut <- qnorm( singlelevel[1] / 2, lower.tail = FALSE)

  sc$lower <- exp( sc$center - zcut * sc$scale )
  sc$upper <- exp( sc$center + zcut * sc$scale )
  if ( npickgene < 0 ) {
    pickedzscore <- abs(z) > zcut
    if( plotit ) {
      rx <- range( sc$x )
      ry <- range( sc$y )
      if( length( negative )) {
        ry <- range( c( ry, exp( y[negative] )))
        if( !meanrank )
          rx <- range( c( rx, exp( intensity[negative] )))
      }
      plot( rx, ry, type = "n", log = logxy, xlab=xlab, ylab=ylab )
      if( main != "" )
        title( main )

      ## critical lines
      lines( sc$x, sc$lower, col = "red" )
      lines( sc$x, sc$upper, col = "red" )
      if( length( singlelevel ) > 1 )
        for( i in seq( 2, length( singlelevel ))) {
          zcut <- qnorm( singlelevel[i] / 2, lower.tail = FALSE)
          lines( sc$x, exp( sc$center - zcut * sc$scale ), lty = 2 )
          lines( sc$x, exp( sc$center + zcut * sc$scale ), lty = 2 )
        }
      ## negligible contrasts
      points( sc$x[!pickedzscore], sc$y[!pickedzscore], cex = 0.25 )
      ## signficant contrasts
      points( sc$x[pickedzscore], sc$y[pickedzscore], cex = 0.75, col = col )
      ## negative values
      if( length( negative ))
        points( exp( intensity[negative] ), exp( y[negative] ),
               cex = 0.25, col = "red" )

      # center line
      lines( sc$x, exp( sc$center ), col = "red", lty = 2 )
    }
    ## get the order and zscore for the top picks
    pickedorder <- seq( n )[pickedzscore]
    z <- z[pickedorder]
    tmp <- order( (-1) * abs( z ))
    pickedorder <- pickedorder[tmp]
    pickedzscore <- z[tmp]
  }
  else {
    pickedorder <- order( (-1) * abs( z ))[1:npickgene]
    pickedzscore <- z[pickedorder]
  }
  ## data frame of picked genes
  fold <- sc$y[pickedorder]
  fold[ fold < 1 ] <- - 1 / fold[ fold < 1 ]
  pick <- data.frame( probe = geneid[pickedorder],
                     average = round( sc$x[pickedorder], 3 ),
                     fold = round( fold, 2 ),
                     pvalue = round( 1 - ( 1 - 2 * pnorm( abs( pickedzscore ),
                       lower.tail = FALSE )) ^ n, 4 ),
                     row.names = NULL )
  ## score structure with centered values
  sc <- as.data.frame( sc )
  sc$y <- ( log( sc$y ) - sc$center ) / sc$scale
  sc$center <- sc$scale <- NULL
  names( sc ) <- c("intensity","score","lower","upper")
  sc$probe <- geneid
  list( pick = pick, score = sc )
}

###########################################################################
## the following function picks genes in the situation of a linear sequence
## of conditions. The argument condi can be numerical (quantitative description
## of the conditions) or ordinal (1 : numberofconditions), depending on the
## situation. The function picks genes according to linear trend (d=1), or
## linear trend and quadratic trend (d=2). Cubic trend can be included if d is
## set to be 3. The argument x is a numberofgenes by numberofconditions matrix
## consisting of logmeasurements.
pickgene.poly <-function( x, # data matrix
                         condi = 1:min(ncol(x),2), # condition levels
                         geneID = NULL,
                         overalllevel = 0.05,
                         npickgene= -1,
                         d=2, # polynomial order (1-3)
                         ylabs = paste( contrastnames, "Trend" ),
                         contrastnames = c("Linear","Quadratic","Cubic"),
                         ... ) {
  if ( d > 3 | d < 1 )
    return( cat( "d should be 1, 2, or 3.\n" ))
  x <- as.matrix( x )
  if( ncol( x ) < 2 )
    return( cat( "x must have at least 2 columns.\n" ))
  d <- min( d, ncol(x)-1 )
  if( is.null( geneID ) | length( geneID ) != nrow( x ))
    geneID <- 1:nrow(x)
  ## get orthogonal polynomial contrasts
  y <- x %*% poly( condi, d )
  intensity <- x %*% rep( 1, length( condi )) / length( condi )

  singlelevel <- adjustlevel( nrow( x ) * d , overalllevel )
  npickgene <- floor( npickgene / d )
  pickedgene <- list()
  pickedgene$pick <- pickedgene$score <- list()

  par( mfrow = c(d,1) )
  for( i in 1:d) {
    tmp <- pickgene.two( y[,i], intensity,
                        geneid = geneID,
                        singlelevel = singlelevel,
                        npick = npickgene,
                        ylab = ylabs[i],
                        ... )
    pickedgene$pick[[i]] <- tmp$pick
    pickedgene$score[[i]] <- tmp$score
  }
  names( pickedgene$pick ) <- contrastnames[1:d]
  names( pickedgene$score ) <- contrastnames[1:d]
  return( pickedgene )
}

###########################################################################
## the following function picks genes in the two-factor (ANOVA) situation.
## It ignores interaction, and picks genes according to trends in each factor.
## Similar functions can easily be written in the same vein for multi-factor
## situation and taking (some) interactions into consideration.
model.pickgene <- function( faclevel,
                           facnames = letters[ seq( length( faclevel )) ],
                           contrasts.fac = "contr.poly",
                           collapse = "+",
                           show = NULL,
                           renorm = 1,
                           modelexpr = formula( paste( "~",
                             paste( facnames, collapse = collapse ))),
                           contrasts.list = contr.list ) {
  ## set up polynomial contrasts
  dd <- expand.grid( apply( as.matrix( rev( faclevel )), 1, seq ))
  names( dd ) <- rev( facnames )
  for( i in facnames )
    dd[[i]] <- factor( dd[[i]] )

  contr.list <- as.list( array( contrasts.fac, length( faclevel )))
  names( contr.list ) <- facnames
  mat <- model.matrix( modelexpr, dd, contrasts = contr.list )
  if( any( renorm != 1 )) {
    if( is.null( show ))
      show <- seq( ncol( mat ) - 1 )
    for( i in seq( length( show )))
      mat[,1+i] <- mat[,1+i] * renorm[show[i]]
  }
  mat
}

###########################################################################
pickgene <- function( data,
                     geneID = 1:nrow(data),
                     overalllevel=0.05,
                     npickgene= -1,
                     marginal = FALSE,
                     rankbased = TRUE, allrank = FALSE,
                     meanrank = FALSE,
                     offset = 0,
                     modelmatrix = model.pickgene( faclevel, facnames,
                       contrasts.fac, collapse, show, renorm ),
                     faclevel = ncol( data ),
                     facnames = letters[ seq( length( faclevel )) ],
                     contrasts.fac = "contr.poly",
                     show = NULL,
                     main = "", renorm = 1, drop.negative = FALSE,
                     plotit = npickgene < 1,
                     mfrow = c(nr,nc), mfcol = NULL,
                     ylab = paste( shownames, "Trend" ),
                     ... ){
  data <- as.matrix( data )
  if( rankbased | allrank ) {
    tmpfn <- function( data ) qnorm( rank( data ) / ( 1 + length( data )))
    data <- if( allrank )
      matrix( tmpfn( data ), nrow( data ), ncol( data ),
             dimnames = dimnames( data ))
    else
      apply( data, 2, tmpfn)
    negative <- numeric(0)
  }
  else { # log transformation, and drop zeroes
    xmin <- apply( data, 1, min, na.rm = TRUE )
    negative <- seq( length( xmin ))[ xmin <= -offset ]
    data <- offset + data
    if( drop.negative ) {
      data <- data[-negative,]
      warn <-  "probes dropped with values below"
    }
    else {
      tmp <- c( data ) > 0
      data[!tmp] <- min( data[tmp] ) / 2
      warn <-  "probes truncated to"
    }
    if( length( negative ))
      warning( paste( length( negative ), warn, offset, "\n" ))
    if( drop.negative )
      negative <- numeric( 0 )
    data <- log( data )
  }

  collapse <- if( marginal ) "+" else "*"
  numcontr <- ncol( modelmatrix ) - 1
  y <- data %*% modelmatrix
  if( is.null( show ))
    show <- seq( numcontr )

  if( meanrank )
    intensity <- rank( y[,1] )
  else
    intensity <- y[,1] / prod( faclevel )
  singlelevel <- adjustlevel( nrow( data ) * numcontr, overalllevel )
  npickgene <- floor ( npickgene / numcontr )

  if( plotit ) {
    nc <- min( 2, ceiling( length( show ) / 3 ))
    nr <- min( 3, ceiling( length( show ) / nc ))
    if( missing( mfcol )) {
      if( !is.null( mfrow ))
        par( mfrow = mfrow )
    }
    else {
      if( !is.null( mfcol ))
        par( mfcol = mfcol )
    }
  }
  pickedgene <- list()
  pickedgene$pick <- pickedgene$score <- list()
  shownames <- dimnames( y )[[2]][1+show]
  main <- array( main, length( shownames ), dimnames = list( shownames ))
  names( ylab ) <- shownames
  for( i in shownames ) {
    tmp <- pickgene.two( y[,i], intensity, geneID, singlelevel,
                        npickgene, main = main[i], plotit = plotit,
                        meanrank = meanrank,
                        negative = negative, ylab = ylab[i], ... )
    pickedgene$pick[[i]] <- tmp$pick
    pickedgene$score[[i]] <- tmp$score
  }
  invisible(pickedgene)
}
