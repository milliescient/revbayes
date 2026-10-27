## name
mvMorphospeciationEvents
## title
Gibbs sampling of the hidden speciation events on a branch
## description
Draws the count of hidden speciation events of one kind on one branch from its full conditional distribution, under dnPhyloMorphospeciation with epsilon.
## details
With epsilon, each branch carries latent Poisson counts of hidden budding and anagenetic events. This is the only move that changes them. Each proposal picks a branch and one kind of event on it and draws the count from its full conditional, so every proposal is accepted.
## authors
June Walker
## see_also
dnPhyloMorphospeciation
mvMorphospeciationAllocation
## example
    epsilon ~ dnBeta(1, 9)
    morph ~ dnPhyloMorphospeciation(tree=bd, epsilon=epsilon)
    morph.clamp(matrix)

    moves.append( mvMorphospeciationEvents(morph, weight=2 * taxa.size()) )
    moves.append( mvMorphospeciationAllocation(morph, weight=taxa.size()) )
    moves.append( mvSlide(epsilon, weight=5) )
## references
