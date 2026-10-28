## name
dnPhyloMorphospeciation
## title
Morphological characters whose changes are speciation events
## description
Discrete morphological characters on a tree drawn from dnFBDSP, under a model in which every change of morphology is a speciation. By default every speciation changes the new species at exactly one character; with `epsilon` it changes each character independently.
## details
Species are morphospecies: a species is a stretch of a lineage over which its morphology does not change. At each speciation, budding or anagenetic, one character chosen uniformly changes to a state drawn uniformly from the others (the Mk kernel). Along a species' stem, from its origination to its first appearance, a character evolves as a Markov chain whose length is set by the rate of hidden speciation, `lambda * p(t) + lambda_a`, divided by the number of characters. The morphology does not change over a species' sampled range. With the characters marginalized, the tree density is that of dnFBDSP.

The number of characters counts every scored character, invariant ones included. For a matrix with its invariant characters dropped, use `coding="variable"`.

Each character has its own number of states K, which sets its change kernel, its starting frequencies 1/K and, under `coding="variable"`, its chance to vary. Left out, `nStates` gives each character as many states as its largest observed state, as `setNumStatesPartition` partitions a matrix for dnPhyloCTMC. One number shares K across every character, which is what an unpartitioned Mk analysis assumes; it takes the number of characters from the data, so it cannot be simulated from before data are clamped. A vector states K for each character, for instance a matrix's declared states where some are never observed.

Which character changed at each speciation is a latent index, sampled by mvMorphospeciationAllocation and reported by `getChangedCharacters()` (0 for the species at the origin). Symmetric speciation is not part of the model.

With `epsilon`, a speciation changes each character independently with probability `epsilon`, so it can change several characters or none. Each branch then carries latent counts of hidden budding and anagenetic events, sampled by mvMorphospeciationEvents. `getNumberOfChanges()`, `getHiddenEventCounts()` and `getHiddenEventMeans()` report them, by node index.

A subtree in which no taxon is coded for any character contributes nothing to the likelihood, so it is never computed, and the hidden counts on its branches are drawn from their Poisson prior without one. Above a taxon coded for no character the counts are summed out, so they are not sampled at all: there `getHiddenEventCounts()` draws them from their Poisson distribution, or with `mean=TRUE` gives its mean.

With `epsilon_a`, an anagenetic speciation changes each character with probability `epsilon_a` and a budding with `epsilon`, so `epsilon / epsilon_a` measures the change carried by budding relative to anagenesis. `getChangeProbabilities()` gives the two probabilities, and `getSpeciationModes()` gives each species' origin: 1 for budding, 2 for anagenesis, 0 at the origin.

With `coding="variable"`, the matrix is taken to hold only the characters that vary among the taxa scored for them, as most published matrices do, and each character's likelihood is divided by its probability of varying. A character is scored on its own taxa, so a taxon missing it does not count toward whether it varies. It needs `epsilon`, since with exactly one change per speciation the characters share their speciations, and dropping the invariant ones does not act on each character alone.
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

    # a matrix that kept only its variable characters
    morph ~ dnPhyloMorphospeciation(tree=bd, epsilon=epsilon, coding="variable")
    morph.clamp(matrix)
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
