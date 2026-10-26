## name
dnPhyloMorphospeciation
## title
Morphological characters whose changes are speciation events
## description
Discrete morphological characters on a tree drawn from dnFBDSP, under a model in which every change of morphology is a speciation. By default every speciation changes the new species at exactly one character; with `epsilon` it changes each character independently.
## details
Species are morphospecies here: a species is a stretch of a lineage over which its morphology does not change. At each speciation, budding or anagenetic, the character that changes is chosen uniformly from all the characters, and its new state is drawn uniformly from the others (the Mk kernel). The parent continues unchanged. Conditional on the tree, a character evolves along a species' stem, from its origination down to its first appearance, as a Markov chain whose length is set by the rates of hidden speciation, `lambda * p(t) + lambda_a`, divided by the number of characters. Over a species' sampled range its morphology does not change, and at the root it starts from the stationary distribution. With the characters marginalized, the tree density is exactly that of dnFBDSP, so the characters redistribute it among morphologies without changing it.

The number of characters counts every scored character, invariant ones included, because a speciation chooses among all of them. Drop invariant characters only if you know how many there were.

Which character changed at each speciation is a latent index per species, sampled by mvMorphospeciationAllocation and reported by `getChangedCharacters()`, which gives 0 for the species that begins at the origin. The rates are read from the dnFBDSP tree, so the characters are recomputed whenever any of them changes. Symmetric speciation is not part of the model.

With `epsilon`, a speciation changes each character independently with probability `epsilon`, so it can change several characters or none. A speciation that changes none is one whose new species differs only in characters that were not scored. The characters are then independent only given how many hidden speciation events each branch holds, so each branch carries a latent count of hidden buddings and one of hidden anagenetic events, Poisson with means `int lambda p dt` and `int lambda_a dt`, sampled by mvMorphospeciationEvents. Given the counts, a character is carried along the branch by `((1 - e) I + e J)^N` for each kind of event. The characters that changed at each observed speciation are sampled by mvMorphospeciationAllocation and their number reported by `getNumberOfChanges()`. `getHiddenEventCounts()` and `getHiddenEventMeans()` give both kinds of event together on the branch above each node, by node index.

With `epsilon_a`, an anagenetic speciation changes each character with probability `epsilon_a`, and a budding with `epsilon`. Expected change per event is `n epsilon` for a budding and `n epsilon_a` for an anagenesis, so `epsilon / epsilon_a` is how much change a budding carries relative to an anagenesis, which is the contrast punctuated equilibrium makes. `getChangeProbabilities()` gives the two probabilities. `getSpeciationModes()` gives each species' origin, 1 for budding and 2 for anagenesis, and 0 for the species at the origin.
## authors
June Walker
## see_also
dnFBDSP
mvMorphospeciationAllocation
mvMorphospeciationEvents
dnPhyloCTMC
## example
    bd ~ dnFBDSP(taxa=taxa, lambda=lambda, mu=mu, psi=psi, rho=rho, lambda_a=lambda_a, origin=origin)
    rec ~ dnFossilRecord(ranges=bd, reporting="complete")
    rec.clamp(taxa)

    morph ~ dnPhyloMorphospeciation(tree=bd)
    morph.clamp(matrix)

    moves.append( mvMorphospeciationAllocation(morph, weight=taxa.size()) )
    changed := morph.getChangedCharacters()

    # any number of characters per speciation, and a separate rate of change at anagenesis
    epsilon ~ dnBeta(1, 9)
    epsilon_a ~ dnBeta(1, 9)
    morph ~ dnPhyloMorphospeciation(tree=bd, epsilon=epsilon, epsilon_a=epsilon_a)
    morph.clamp(matrix)

    moves.append( mvMorphospeciationEvents(morph, weight=2 * taxa.size()) )
    moves.append( mvMorphospeciationAllocation(morph, weight=taxa.size()) )
    moves.append( mvSlide(epsilon, weight=5) )
    moves.append( mvSlide(epsilon_a, weight=5) )
    ratio := epsilon / epsilon_a
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
