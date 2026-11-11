## name
mvExtinctionRJSwitch
## title
Jump a species' extinction time to and from the present
## description
This move switches one species' extinction time between the present, representing unsampled survival, and an age older than the present.
## details
When `rho` < 1, a species reported extinct may have survived without being sampled. Its extinction time therefore has a point mass at the present, which continuous moves never propose. This move switches between that point mass and an older extinction age. Species reported extant are never chosen, and the move does nothing when `rho` = 1.

You can use this move with the range matrix of a `dnFBDRP` or `dnBDS`, or with the tree of a `dnFBDSP`.
## authors
June Walker
## see_also
mvSingletonRJSwitch
mvMatrixElementSlide
mvRJSwitch
## example
    bd ~ dnFBDRP(taxa=taxa, lambda=lambda, mu=mu, psi=psi, rho=0.5)

    moves.append( mvMatrixElementSlide(bd, delta=1, weight=taxa.size()) )
    moves.append( mvExtinctionRJSwitch(bd, weight=taxa.size()/5) )
## references
