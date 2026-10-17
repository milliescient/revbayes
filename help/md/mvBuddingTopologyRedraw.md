## name
mvBuddingTopologyRedraw

## title
Gibbs draw of a budding topology

## description
Redraws the whole budding (asymmetric speciation) topology of a `dnFBDSP` tree, holding the species ranges fixed.

## details
Each lineage is assigned an ancestor drawn uniformly from those alive at its birth. Given the ranges, every compatible tree has the same probability, so the draw is from the exact conditional and every proposal is accepted. With character data on the tree the move is still valid, but acceptance then rests on the likelihood, so pair it with `mvFNPR` or `mvNNI`.

The move requires pure budding and refuses to be set up when `lambda_a > 0`.

It never changes which species buds at a given birth age, so pair it with `mvRotateNode`. Ranges and node ages are untouched: use `mvFossilTipTimeUniform` for the extinction times and `mvNodeTimeSlideUniform` with `mvRootTimeSlideUniform` for the speciation times.

## authors
June Walker

## see_also
dnFossilizedBirthDeathSpeciation
mvRotateNode
mvFNPR

## example
tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, timeline=timeline, taxa=taxa)
rec ~ dnFossilRecord(ranges=tr, reporting="firstlast")
rec.clamp(taxa)

moves.append( mvBuddingTopologyRedraw(tr, weight=taxa.size()) )
# relabels which lineage carries each birth age; the Gibbs move cannot
moves.append( mvRotateNode(tr, weight=taxa.size()) )
moves.append( mvFossilTipTimeUniform(tr, weight=taxa.size()) )
moves.append( mvNodeTimeSlideUniform(tr, weight=taxa.size()) )
moves.append( mvRootTimeSlideUniform(tr, origin, weight=2) )

## references
