\name{pickgene-internal}
\alias{pickgene-internal}
\alias{adjustlevel}
\alias{chen.poly}
\alias{chipnorm}           
\alias{dencont}
\alias{dencum}
\alias{denlines}           
\alias{do.oddsplot}
\alias{fitgg}
\alias{gammaden}           
\alias{holms}
\alias{lod.ggb}
\alias{lod.plot}           
\alias{lodprobes}
\alias{loglik}
\alias{makecont}           
\alias{multipickgene}
\alias{nlminb}
\alias{nloglik}            
\alias{normal.richmond}
\alias{npdiag}
\alias{nploglik}           
\alias{orangene}
\alias{pickedchisq}
\alias{pickedhist}         
\alias{pickedpair}
\alias{pickedscore}
\alias{pickgene2}          
\alias{pickgene.poly}
\alias{pickgene.two}
\alias{pmarg}              
\alias{predrecur}
\alias{rangene}
\alias{rankgene}           
\alias{robustbox}
\alias{s.check0}
\alias{s.check1}           
\alias{s.check2}
\alias{shrinkplot}
\alias{sixden}             
\alias{s.marg}
\alias{toprankgene}
\alias{twoarray.norm}      
\alias{twoarray.plot}
\alias{twowayanovapickgene}
\title{Internal pickgene functions.}
\description{  These are generally not to be called by the user.
}
\usage{
adjustlevel(ntest, alpha)
chen.poly(cv, err)
chipnorm(xx, chip)
dencont(x, y, align, crit, xlim, ylim, dolog, byranks, dif,
                 ave, numlines, levels.z)
dencum(x, y, align, crit, xlim, ylim, dolog, byranks,
                 standardize, dif, ave, splineit, numlines, show,
                 levels.z)
denlines(x, y, align, crit, xlim, ylim, dolog, dif, ave,
                 numlines, offset)
do.oddsplot(data, main, theta, col, redo, conditions, identifier,
                 ...)
fitgg(xx, yy, start)
gammaden(x, a, b)
holms(x, alpha, cut)
lod.ggb(x, y, theta)
lod.plot(data, x, y, theta, filename, probe, xlab, ylab, ps,
                 col, lowlod, ...)
lodprobes(xx, yy, theta, lod, probes, col, lowlod, offset)
loglik(theta, xx, yy)
makecont(x, y, size, cex, levels)
multipickgene(...)
nlminb(start, objective, lower, xx, yy, zz, use.optim)
nloglik(theta, xx, yy)
normal.richmond(foo, channel)
npdiag(xx, yy, aa, a0, nu, pp)
nploglik(theta, xx, yy, zz)
orangene(n, center, spread, contamination, alpha, noise,
                 omega)
pickedchisq(pick, show, title, plotit, alpha)
pickedhist(pick, show, title, p1, plotit, rotate, mfrow, bw)
pickedpair(x, columns, description, probe, renorm, pick, main,
                 ...)
pickedscore(pick, description, show, alpha, xlab, ylab, main,
                 mfrow)
pickgene2(...)
pickgene.poly(x, condi, geneID, overalllevel, npickgene, d, ylabs,
                 contrastnames, ...)
pickgene.two(y, intensity, geneid, singlelevel, npickgene,
                 meanrank, xlab, ylab, main, plotit, col, negative,
                 ...)
pmarg(xx, yy, theta, nsupp)
predrecur(xx, theta, gridlim)
rangene(n, center, spread, contamination, alpha, noise,
                 omega)
rankgene(xx, yy, fits)
robustbox(y, x, nslice, xlab, ylab, shrink, crit, overalllevel,
                 cex, lwd, plotit)
s.check0(xx, yy, theta1, theta2, chip)
s.check1(xx, yy, theta, chip)
s.check2(foo, xa, ya, thetaa, xb, yb, thetab, spots)
shrinkplot(xx, yy, fits, chip)
sixden(x, y, align, crit, xlim, dolog, dif, ave)
s.marg(xx, yy, aa, a0, nuA, nu0, p)
toprankgene(yy, n)
twoarray.norm(foo, ..., conditions, reduce, identifier)
twoarray.plot(mydata, main, theta, conditions, identifier)
twowayanovapickgene(x, fac1level, fac2level, ...)
}
\author{ Brian S. Yandell, yandell@stat.wisc.edu }
\keyword{internal}
