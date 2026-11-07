## name
dnFossilizedBirthDeathSpeciation
## title
The fossilized birth-death range process over trees
## description
The fossilized birth death speciation process (FBDSP) describes the diversification and sampling of extant and extinct species trees under a mixed model of asymmetric (budding) and anagenetic speciation.
## details
Fossil species are represented by a collection of fossil occurrences with uncertainty. Speciation, extinction and sampling rates may be time-homogeneous or piecewise time-heterogeneous. If time-heterogeneous rates are provided, then a vector of rate change time-points must also be provided. Like `dnFBDRP`, this is the birth-death range process on its own (here over trees); pair it with a `dnFossilRecord` node clamped to the occurrences to add the probability of the fossil record. The first and last appearance of each species are latent and are not carried by the tree, so sample them with `mvStratigraphicRange`. By default the process runs to age zero. `present` moves that boundary, so a record truncated at a stratigraphic horizon is analysed in the ages it already carries: the timeline starts there, no occurrence may be younger, and `rho` applies at that age rather than at zero.

The tree is an extended tree: each extinct species' tip is its extinction time, below its last appearance, and an extant species' tip is at the present. Every age in the tree needs a move. Tips are sampled by `mvFossilTipTimeUniform` between the present and the youngest occurrence. `mvNodeTimeSlideUniform` skips the root, so pair it with `mvRootTimeSlideUniform`.

Fossil sampling rates can vary among species, as in `dnFBDRP`. Each species draws its own component, whether it arose by budding, symmetric speciation or anagenesis. `mvPsiMixtureAllocation` samples the components and `getPsiMixtureComponents()` reports them.

With `beta > 0` a speciation may also be symmetric (bifurcating): the ancestral species ends and both daughters begin new ones.

A sampled ancestor's successor arose by anagenesis or by a symmetric speciation whose other daughter left no sample. The successor's continuation flag holds which: on where its lineage did not bifurcate away, so the speciation mode is read from the flags at every node. A Newick carries it as `mode=A` (anagenesis) or `mode=S` (symmetric) on the successor, and a successor without it is taken as anagenetic where lambda_a > 0.

`getSampledTree()` returns the sampled tree the extended tree projects to: each extinct species cut back to its last appearance, with a sampled ancestor there where its lineage leads on, and each node's species label as `sp`. Its tips are last appearances, so they move with `mvStratigraphicRange` and not with the tip moves. Log it with `mnFile` for tree summaries; it is computed only when something reads it, so logging it does not slow the chain.

`getObservedSpeciation(mode)` and `getHiddenSpeciation(mode)` return per-branch counts, indexed as dnPhyloCTMC's `branchRates` are. The observed count is the event that began the branch's species; the hidden count is the expected number of unobserved speciations along the branch. `mode` is `"budding"`, `"bifurcation"`, `"anagenesis"`, `"cladogenesis"` (budding and bifurcation), or `"all"` (the default).

## authors
June Walker
## see_also
dnFossilizedBirthDeathRange
dnFossilRecord
mvPsiMixtureAllocation
mvFossilTipTimeUniform
mvRootTimeSlideUniform
mvStratigraphicRange
## example
lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)

min_age = 0.0
for(i in 1:taxa.size())
{
	if ( taxa[i].getMinAge() > min_age )
	{
		min_age = taxa[i].getMinAge()
	}
}

origin_offset ~ dnExp(1/10)
origin := min_age + origin_offset
moves.append( mvSlide(origin_offset, weight = 2) )

bd ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

# the fossil record, conditioned on the range process and clamped to the observed occurrences
rec ~ dnFossilRecord(ranges=bd, reporting="firstlast")
rec.clamp(taxa)

moves.append( mvFNPR(bd, weight = taxa.size()) )
moves.append( mvNodeTimeSlideUniform(bd, weight = taxa.size()) )
moves.append( mvRootTimeSlideUniform(bd, origin=origin, weight = taxa.size()) )
moves.append( mvFossilTipTimeUniform(bd, weight = taxa.size()) )

# which lineage continues its ancestor's species at each speciation event
moves.append( mvSpeciationMode(bd, weight = taxa.size()) )

# the first and last appearances, which the tree does not carry
moves.append( mvStratigraphicRange(bd, weight = taxa.size()) )

# the sampled tree, tips at last appearances, for tree summaries
sampled := bd.getSampledTree()
monitors.append( mnFile(sampled, filename="output/sampled.trees", printgen=100) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
