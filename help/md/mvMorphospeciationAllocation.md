## name
mvMorphospeciationAllocation
## title
Gibbs sampling of the character that changed at a speciation
## description
Draws the character that changed at one species' speciation from its full conditional distribution.
## details
Under dnPhyloMorphospeciation, each species carries the index of the character that changed when it arose. The index is not part of the value, so this is the only move that changes it. Each proposal picks a species with a speciation event at random and draws its index given the tree and every species' morphology, from two passes over the tree, so every proposal is accepted. With `epsilon`, a species carries a set of changed characters instead, and the same two passes give each character's odds of being in it. Each character is then drawn on its own.
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
