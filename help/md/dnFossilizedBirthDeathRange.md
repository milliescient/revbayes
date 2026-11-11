## name
dnFossilizedBirthDeathRange
## title
The fossilized birth-death range process
## description
This distribution describes a matrix of species ranges under budding speciation, extinction and fossil sampling, allowing for species that left no samples.
## details
Each row contains one species' extinction time, last appearance, first appearance and birth time, in that order. The first and last appearances lie within the reported occurrence ages. Use matrix element moves to sample all four times. Use `mvMatrixElementSwap` with `col=4` to exchange birth times between species and explore birth orderings.

Add a `dnFossilRecord` node with `ranges` set to this matrix and clamp it to the taxa. This node includes the fossil occurrences and sets the reporting model.

`lambda`, `mu` and `psi` take one value or one per interval defined by `timeline`. Without `origin`, the process starts at the oldest birth. `present` defaults to 0. Set it to analyse a record cut off at a stratigraphic horizon. Choose `condition` to match how the data were selected.

Each birth is weighted by the number of coexisting sampled species, and a birth with no such species has zero probability.

RevBayes warns at the start of an MCMC about missing required moves, including `mvExtinctionRJSwitch` when `rho` < 1.

`getOriginationTimes()` draws each species' own origination time. This differs from its birth time when it descends from an unsampled lineage.
## authors
June Walker
## see_also
dnFossilRecord
mvMatrixElementSwap
mvExtinctionRJSwitch
mvSingletonRJSwitch
## example
    taxa <- readTaxonData("fossil_ranges.tsv")

    lambda ~ dnExponential(10)
    mu     ~ dnExponential(10)
    psi    ~ dnExponential(10)

    bd ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

    rec ~ dnFossilRecord(ranges=bd, reporting="firstlast")
    rec.clamp(taxa)

    n <- taxa.size()
    moves.append( mvMatrixElementScale(bd, weight=n) )
    moves.append( mvMatrixElementSlide(bd, weight=n) )
    moves.append( mvMatrixElementSwap(bd, col=4, weight=n/5) )

## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
	- citation: Improved estimation of macroevolutionary rates from fossil data using a Bayesian framework. Silvestro, Daniele et al. Paleobiology, 45:546-570.
	  doi: https://doi.org/10.1017/pab.2019.23
	  url: https://www.cambridge.org/core/journals/paleobiology/article/improved-estimation-of-macroevolutionary-rates-from-fossil-data-using-a-bayesian-framework/334F08A74A6C92F1FEAD91A71FE59A1C
