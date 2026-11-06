## name
mvExtinctionTimeGibbs

## title
Gibbs move for extinction times

## description
On an extended dnFBDSP tree, redraws one fossil tip's extinction time from its distribution given the rest of the model.

## details
On an extended tree each extinct tip is its species' extinction. Given the tree, the rates and the species' last appearance, that time has a conditional density proportional to q~(tau_K -> d) times the probability of the range's end at d, between the present and the younger of the last appearance and the tip's parent. The move draws from it exactly, so with one fossil sampling rate it is always accepted. Under a psi mixture the draw is an independence proposal, still correct, accepted at a lower rate. It moves one tip per proposal and skips sampled ancestors, which sit at their nodes. On a sampled tree (extended=FALSE) the extinction times are integrated out, and the move refuses to be built.

## authors
June Walker

## see_also
dnFossilizedBirthDeathSpeciation
mvFossilTipTimeUniform

## example
tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa, extended=TRUE)
moves.append( mvExtinctionTimeGibbs(tr, weight=taxa.size()/5) )

## references
