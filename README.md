clhs.jl
=======

Implementation of conditioned Latin Hypercube Sampling in Julia
---------------------------------------------------------------

This is a port of [my R package `clhs`](http://cran.r-project.org/web/packages/clhs/index.html) to Julia.

It implements the conditioned Latin Hypercube Sampling, as [published by Minasny and McBratney (2006)](http://www.sciencedirect.com/science/article/pii/S009830040500292X). This method proposes to stratify sampling in presence of ancillary data. 
