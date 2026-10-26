## name
mvMorphospeciationEvents
## title
Gibbs sampling of the hidden speciation events on a branch
## description
Draws the count of hidden speciation events of one kind on one branch from its full conditional distribution, under dnPhyloMorphospeciation with epsilon.
## details
With epsilon, a speciation can change any number of characters, and the characters are independent only given how many speciation events each branch holds. Each branch therefore carries a latent count of hidden buddings and one of hidden anagenetic events, each Poisson with its integrated rate. The counts are not part of the value, so this is the only move that changes them. Each proposal picks a branch that can hold events and one kind of event on it, and draws that count given the other, the tree and every species' morphology. Every character's likelihood is linear in the branch's transition matrix, so two passes over the tree give the full conditional for every count, and every proposal is accepted.
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
