## name
mvSpeciationMode

## title
Speciation mode move

## description
Moves the continuation of an ancestral species between the children of a random internal node, leaving every node age and the topology unchanged. On a sampled tree whose per-branch event counts are read, it also redraws where a sampled ancestor's successor changed species. Also available as mvRotateNode.

## details
Under budding (asymmetric) speciation exactly one child of each node continues its ancestor's species, and the others begin new species there. Which one continues is a free parameter, separate from the topology and the node ages, and this move samples it. The new continuing child is drawn uniformly among the others.

Where beta > 0, a node may also have no continuing child: its species ends there and both children begin new ones. The move then draws uniformly among the node's other states, one per child plus symmetric, so the proposal is its own reverse.

At a sampled ancestor on an extended tree the tip carries the ancestor's species, which ends there, and the other child is its successor. The successor's flag says whether its lineage bifurcated away from the ancestor's: off for a symmetric speciation whose other daughter left no sample, on for an anagenesis, where the lineage carries on under a new species. Where both lambda_a and beta are positive, each successor joins the nodes the move picks from, and picking one flips its mode. Otherwise the mode is fixed by the rate that is positive.

On a sampled tree (extended=FALSE), a sampled ancestor's named successor first appears below it, so its lineage changed species at least once on the way down, unobserved. Once any of the per-branch event accessors (getHiddenBuddings and the others) is created, dnFBDSP holds the first of those changes as a latent: its age, and whether it was an anagenesis, a budding whose parent side left no sample, or a symmetric speciation whose other daughter left none. The accessors report it as that branch's observed event, and count hidden events only below it, so the same model reads both tree types. Each named successor then joins the nodes the move picks from, and picking one redraws its first change exactly from its density given the tree, so only what reads the counts decides its acceptance. Without the accessors there is no latent and the move only moves continuations.

## authors
June Walker

## see_also
dnFossilizedBirthDeathSpeciation
mvBuddingTopologyRedraw
mvStratigraphicRange

## example
tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, timeline=timeline, taxa=taxa)

# which lineage continues the ancestral species at each speciation event
moves.append( mvSpeciationMode(tr, weight=taxa.size()/2) )

# with beta > 0 it also moves nodes between budding and symmetric
beta ~ dnBeta(1, 1)
tr_sym ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, beta=beta, taxa=taxa)
moves.append( mvSpeciationMode(tr_sym, weight=taxa.size()/2) )

## references
