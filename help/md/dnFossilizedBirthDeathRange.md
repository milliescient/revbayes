## name
dnFossilizedBirthDeathRange
## title
The fossilized birth-death range process
## description
The fossilized birth death range process (FBDRP) describes the distribution of a matrix of species divergence and extinction times under a model of asymmetric speciation and sampling of extinct species.
## details
This distribution is the birth-death range process on its own. Pair it with a `dnFossilRecord` node clamped to the occurrences to add the probability of the fossil record; that node also selects the reporting model. Without one, the process does not condition on the fossil record.

The `condition` argument has to match whatever the data were selected for. `time` conditions on nothing, `sampling` on at least one sampled taxon, `survival` on at least one lineage alive at the present, and `nTaxa` on the number of taxa the record reports as extant, which is what data simulated under a fixed number of extant species needs.

Fossil species are represented by a collection of fossil occurrences with uncertainty. Speciation, extinction and sampling rates may be time-homogeneous or piecewise time-heterogeneous. If time-heterogeneous rates are provided, then a vector of rate change time-points must also be provided. Under the hood, the fossil data is augmented with the first and last appearance of each species. These are latent, held in columns 2 and 3 of the value, and are sampled during MCMC by the matrix element moves. By default the process runs to age zero. `present` moves that boundary, so a record truncated at a stratigraphic horizon is analysed in the ages it already carries: the timeline starts there, no occurrence may be younger, and `rho` applies at that age rather than at zero.

`coexistence` weights each birth by how many lineages were alive to produce it, and defaults on. Turning it off matches `dnBDS`.

The process exposes its latent ages as member methods. `getDivergenceTimes` gives the time each range attaches to the nearest sampled range. `getOriginationTimes` draws each species' own origination from its conditional, or with `mean=TRUE` returns the conditional mean. `getExtantStatus` and `getSingletonStatus` are indicators, so their posterior means are probabilities; `getSingletonStatus` is 1 where a taxon's first and last appearances are one age. `getExtinctionTimes` reports the extinction times.

Fossil sampling rates can also vary among species. With `psiMixtureRates`, each species draws a mixture component at birth and samples fossils at that component's rate times `psi`. `psiMixtureWeights` gives the component probabilities and defaults to equal weights. `mvPsiMixtureAllocation` samples the components and `getPsiMixtureComponents()` reports them.

## authors
June Walker
## see_also
dnFossilRecord
mvStratigraphicRange
mvPsiMixtureAllocation
fnDiscretizeGamma
## example
lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)

bd ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

# the fossil record, conditioned on the range process and clamped to the observed occurrences
rec ~ dnFossilRecord(ranges=bd, reporting="firstlast")
rec.clamp(taxa)

moves.append( mvMatrixElementScale(bd, weight=taxa.size()) )
moves.append( mvMatrixElementSlide(bd, weight=taxa.size()) )

# fossil sampling rates that vary among species, as a discretized gamma
alpha ~ dnExp(1)
bd_mix ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa,
                 psiMixtureRates=fnDiscretizeGamma(alpha, alpha, 4))
moves.append( mvPsiMixtureAllocation(bd_mix, weight=taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
	- citation: Improved estimation of macroevolutionary rates from fossil data using a Bayesian framework. Silvestro, Daniele et al. Paleobiology, 45:546-570.
	  doi: https://doi.org/10.1017/pab.2019.23
	  url: https://www.cambridge.org/core/journals/paleobiology/article/improved-estimation-of-macroevolutionary-rates-from-fossil-data-using-a-bayesian-framework/334F08A74A6C92F1FEAD91A71FE59A1C
