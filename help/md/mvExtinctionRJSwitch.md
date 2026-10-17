## name
mvExtinctionRJSwitch
## title
Jump a taxon's extinction time to and from the present
## description
Moves one taxon's extinction time between the present, where it survived without being sampled, and a time above it.
## details
When rho is below one, a taxon reported extinct may have survived unsampled, so its extinction time has a point mass at the present. The continuous element moves never propose that point, so this move is needed to sample it. Taxa reported extant are never chosen.

At rho = 1 the move has nothing to do.
## authors
June Walker
## see_also
mvMatrixElementSlide
mvRJSwitch
## example
    bd ~ dnFBDRP(taxa=taxa, lambda=lambda, mu=mu, psi=psi, rho=0.5, timeline=timeline)

    moves.append( mvMatrixElementSlide(bd, delta=1, weight=taxa.size()/10) )
    moves.append( mvExtinctionRJSwitch(bd, weight=taxa.size()/10) )
## references
