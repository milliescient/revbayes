## name
mvMorphospeciationAllocation
## title
Gibbs sampling of the character that changed at a speciation
## description
Draws the character that changed at one species' speciation from its full conditional distribution.
## details
Under dnPhyloMorphospeciation, each species carries the index of the character that changed when it arose, or with `epsilon` a set of changed characters. This is the only move that changes them. Each proposal picks a speciation at random and draws from the full conditional, so every proposal is accepted.
## authors
June Walker
## see_also
dnPhyloMorphospeciation
mvMorphospeciationEvents
mvGibbsMixtureAllocation
## example
    morph ~ dnPhyloMorphospeciation(tree=bd)
    morph.clamp(matrix)

    moves.append( mvMorphospeciationAllocation(morph, weight=taxa.size()) )
## references
