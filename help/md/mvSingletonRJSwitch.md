## name
mvSingletonRJSwitch
## title
Jump a taxon between one appearance age and two
## description
Expands a taxon's single appearance age into a first and a last appearance within its reported interval, or collapses the two into one.
## details
Under reporting="presence" a taxon with one reported occurrence over an interval may have more, so its first and last appearances may differ within that interval. The continuous moves never change that, so this move is needed to sample it. The last appearance stays where it is: an expansion draws the first appearance uniformly above it, within the interval and below the divergence, and a collapse sets the first appearance onto it. Taxa with more than one occurrence, or an exact age, are never chosen.

The move takes the range matrix of a dnFBDRP or the tree of a dnFBDSP. Under the other reporting models it has nothing to do. `getSingletonStatus()` on the process reports each taxon's state, 1 for a single age, so its posterior mean is the probability that the reported occurrence was the only one.
## authors
June Walker
## see_also
mvExtinctionRJSwitch
dnFossilRecord
## example
    bd ~ dnFBDRP(taxa=taxa, lambda=lambda, mu=mu, psi=psi, rho=1, timeline=timeline)
    rec ~ dnFossilRecord(ranges=bd, reporting="presence")
    rec.clamp(taxa)

    moves.append( mvMatrixElementSlide(bd, delta=1, weight=taxa.size()/10) )
    moves.append( mvSingletonRJSwitch(bd, weight=taxa.size()/10) )
## references
