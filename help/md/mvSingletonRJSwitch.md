## name
mvSingletonRJSwitch
## title
Jump a taxon between one appearance age and two
## description
Expands a taxon's single appearance age into a first and a last appearance within its reported interval, or collapses the two into one.
## details
Under reporting="presence" a taxon with one reported occurrence over an interval may have more, so its first and last appearances may differ. The continuous moves cannot change that. An expansion draws the first appearance uniformly above the last, within the interval and below the divergence; a collapse sets it onto the last. Taxa with more than one occurrence, or an exact age, are never chosen.

The move takes the range matrix of a dnFBDRP or the tree of a dnFBDSP, and does nothing under the other reporting models. `getSingletonStatus()` reports each taxon's state.
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
