## name
dnPhyloMorphospeciation
## title
Morphological characters whose changes are speciation events
## description
Discrete morphological characters on a tree drawn from dnFBDSP, under a model in which every speciation changes the new species at exactly one character, and every change of morphology is a speciation.
## details
Species are morphospecies here: a species is a stretch of a lineage over which its morphology does not change. At each speciation, budding or anagenetic, the character that changes is chosen uniformly from all the characters, and its new state is drawn uniformly from the others (the Mk kernel). The parent continues unchanged. Conditional on the tree, a character evolves along a species' stem, from its origination down to its first appearance, as a Markov chain whose length is set by the rates of hidden speciation, `lambda * p(t) + lambda_a`, divided by the number of characters. Over a species' sampled range its morphology does not change, and at the root it starts from the stationary distribution. With the characters marginalized, the tree density is exactly that of dnFBDSP, so the characters redistribute it among morphologies without changing it.

The number of characters counts every scored character, invariant ones included, because a speciation chooses among all of them. Drop invariant characters only if you know how many there were.

Which character changed at each speciation is a latent index per species, sampled by mvMorphospeciationAllocation and reported by `getChangedCharacters()`, which gives 0 for the species that begins at the origin. The rates are read from the dnFBDSP tree, so the characters are recomputed whenever any of them changes. Symmetric speciation is not part of the model.
## authors
June Walker
## see_also
dnFBDSP
mvMorphospeciationAllocation
dnPhyloCTMC
## example
    bd ~ dnFBDSP(taxa=taxa, lambda=lambda, mu=mu, psi=psi, rho=rho, lambda_a=lambda_a, origin=origin)
    rec ~ dnFossilRecord(ranges=bd, reporting="complete")
    rec.clamp(taxa)

    morph ~ dnPhyloMorphospeciation(tree=bd)
    morph.clamp(matrix)

    moves.append( mvMorphospeciationAllocation(morph, weight=taxa.size()) )
    changed := morph.getChangedCharacters()
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
