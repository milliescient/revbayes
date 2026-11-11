## name
mvBuddingTopologyRedraw
## title
Gibbs draw of a budding topology
## description
This move redraws the whole topology of a pure-budding `dnFBDSP` tree while holding the species ranges fixed. It does not change ages.
## details
The move attaches each species to an ancestor drawn uniformly from the species alive at its birth. Given the ranges, every compatible tree has the same probability, so every proposal is accepted without character data. With character data, use it alongside `mvFNPR`.

The move does not change which daughter continues at a node, so use it with `mvSpeciationMode`. Setting up this move with `lambda_a` > 0 or `beta` > 0 is an error.
## authors
June Walker
## see_also
dnFossilizedBirthDeathSpeciation
mvSpeciationMode
mvFNPR
## example
    tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)
    rec ~ dnFossilRecord(ranges=tr, reporting="firstlast")
    rec.clamp(taxa)

    moves.append( mvBuddingTopologyRedraw(tr, weight=4) )
    moves.append( mvSpeciationMode(tr, weight=taxa.size()) )
    moves.append( mvFossilTipTimeUniform(tr, origin=origin, weight=taxa.size()) )
    moves.append( mvNodeTimeSlideUniform(tr, weight=taxa.size()) )
    moves.append( mvRootTimeSlideUniform(tr, origin=origin, weight=2) )
## references
