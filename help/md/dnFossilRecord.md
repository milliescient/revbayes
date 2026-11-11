## name
dnFossilRecord
## title
The fossil record of a set of stratigraphic ranges
## description
This distribution gives the probability of observed fossil occurrences given species ranges from `dnFBDRP`.
## details
Set `ranges` to the range process and clamp this node to the same taxa. The occurrences, `psi` and the timeline are read from the range process. Without this node, the likelihood of the occurrences is not computed, and RevBayes warns at the start of an MCMC.

Set `reporting` to describe how sampled occurrences entered the data. Use `"complete"` when every sampled occurrence is reported. The default, `"firstlast"`, treats the oldest and youngest occurrence of each species as reported, with an unknown count between them. Use `"presence"` when each reported interval contains at least its reported count. With `"presence"`, intervals for the same species must not overlap.
## authors
June Walker
## see_also
dnFossilizedBirthDeathRange
mvSingletonRJSwitch
## example
    taxa <- readTaxonData("fossil_ranges.tsv")

    lambda ~ dnExponential(10)
    mu     ~ dnExponential(10)
    psi    ~ dnExponential(10)

    bd ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

    rec ~ dnFossilRecord(ranges=bd, reporting="firstlast")
    rec.clamp(taxa)

    moves.append( mvMatrixElementScale(bd, weight=taxa.size()) )
    moves.append( mvMatrixElementSlide(bd, weight=taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
