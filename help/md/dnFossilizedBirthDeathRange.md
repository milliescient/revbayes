## name
dnFossilizedBirthDeathRange
## title
The fossilized birth-death range process
## description
The fossilized birth death range process (FBDRP) describes the distribution of a matrix of species divergence and extinction times under a model of asymmetric speciation and sampling of extinct species.
## details
This distribution is the birth-death range process on its own: it describes how species ranges diversify, and the fossil occurrences bound those ranges, but it does not include the probability of the fossil record itself. Pair it with a `dnFossilRecord` node clamped to the occurrences to get the full model. That node also selects the reporting model, i.e. how sampled occurrences reach the record: `reporting="complete"` (all reported) or `reporting="firstlast"` (only the oldest and youngest). A range process with no `dnFossilRecord` attached is a valid model, but it does not condition on the fossil record at all.

The `condition` argument has to match whatever the data were selected for. `time` conditions on nothing, `sampling` on at least one sampled taxon, `survival` on at least one lineage alive at the present, and `nTaxa` on the number of taxa the record reports as extant, which is what data simulated under a fixed number of extant species needs.

Fossil species are represented by a collection of fossil occurrences with uncertainty. Speciation, extinction and sampling rates may be time-homogeneous or piecewise time-heterogeneous. If time-heterogeneous rates are provided, then a vector of rate change time-points must also be provided. Under the hood, the fossil data is augmented with the first and last appearance of each species. These are latent, held in columns 2 and 3 of the value, and are sampled during MCMC by the matrix element moves. By default the process runs to age zero. `present` moves that boundary, so a record truncated at a stratigraphic horizon is analysed in the ages it already carries: the timeline starts there, no occurrence may be younger, and `rho` applies at that age rather than at zero.

The process exposes its latent ages as member methods. `getDivergenceTimes` gives the time each range attaches to the nearest range the record samples, which is not the same as where the species itself began, because ancestors that left no occurrence were never named and a chain of them may stand in between. `getOriginationTimes` gives the species' own origination, drawn from its exact conditional, so two calls on one state return two different draws; pass `mean=TRUE` to integrate that conditional instead and get its mean without simulation noise. `getExtantStatus` and `getExtinctionTimes` report the state as it stands, so the first is an indicator you average over the posterior to get a probability rather than a probability itself.

Fossil sampling rates can also vary among species. With `psiMixtureRates`, each species draws a mixture component at birth, independently of its parent, and samples fossils at that component's rate times `psi`. `psiMixtureWeights` gives the probability of each component and defaults to equal weights. `psi` is then the mean sampling rate when the rates average one under the weights, as `fnDiscretizeGamma(alpha, alpha, K)` does. Each species carries its component as a latent variable, which `mvPsiMixtureAllocation` samples and `getPsiMixtureComponents()` reports. The probability that a lineage leaves no sampled descendant accounts for the mixture, because unsampled descendants draw their own components.

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
