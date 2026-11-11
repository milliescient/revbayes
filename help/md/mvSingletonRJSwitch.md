## name
mvSingletonRJSwitch
## title
Jump a species between one appearance age and two
## description
This move expands a species' single appearance into a first and last appearance within its reported interval, or collapses the two into one.
## details
Under `reporting="presence"`, a species with one reported occurrence over an interval may have had more. An expansion draws the first appearance uniformly above the last, within the interval and below the birth time. A collapse sets the first appearance onto the last. Species with more than one occurrence, or an exact age, are never chosen.

You can use this move with the range matrix of a `dnFBDRP` or `dnBDS`, or with the tree of a `dnFBDSP`. It does nothing under other reporting models.
## authors
June Walker
## see_also
mvExtinctionRJSwitch
dnFossilRecord
## example
    bd ~ dnFBDRP(taxa=taxa, lambda=lambda, mu=mu, psi=psi, rho=1)
    rec ~ dnFossilRecord(ranges=bd, reporting="presence")
    rec.clamp(taxa)

    moves.append( mvMatrixElementSlide(bd, delta=1, weight=taxa.size()) )
    moves.append( mvSingletonRJSwitch(bd, weight=taxa.size()/5) )
## references
