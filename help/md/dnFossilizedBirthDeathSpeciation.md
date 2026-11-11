## name
dnFossilizedBirthDeathSpeciation

## title
Fossilized birth-death speciation process

## description
This distribution describes a time tree of species under the fossilized birth-death range process. You provide species and their fossil occurrences in `taxa`, with each occurrence specified as an age or an age interval.

## details
You must provide `origin`. Rates can take one value or one value per interval defined by `timeline`. This distribution accounts for the birth-death range process only. To include the fossil occurrences, add a `dnFossilRecord` node with `ranges` set to this tree and clamp it to `taxa`.

The process allows budding, symmetric speciation and anagenesis. In budding, the ancestral species continues and one daughter is new. In symmetric speciation, the ancestral species ends and both daughters are new. In anagenesis, a species becomes a new species along a lineage. Speciation is symmetric with probability `beta` and budding otherwise. Anagenesis occurs at rate `lambda_a`. Both default to 0, giving pure budding.

Each tip represents a species' extinction time, below its last appearance. An extant species' tip is at the present. A fossil species whose lineage continues past it is a sampled ancestor, possible only when `lambda_a` > 0 or `beta` > 0. First and last appearances are not part of the tree. For tree summaries, `getSampledTree()` cuts the tree back to each species' last appearance.

Besides the usual tree moves, you need `mvSpeciationMode` to change which daughter continues the ancestral species and whether a node is budding or symmetric. No other move changes these assignments. You also need `mvStratigraphicRange` for first and last appearances. When `lambda_a` > 0 or `beta` > 0, you also need `mvCollapseExpandFossilBranch` and `mvWilsonBalding`. RevBayes warns at the start of an MCMC when a needed move is missing. The example shows a full set.

## authors
June Walker

## see_also
dnFossilizedBirthDeathRange
dnFossilRecord
mvSpeciationMode
mvStratigraphicRange
mvCollapseExpandFossilBranch
mvWilsonBalding
mvPsiMixtureAllocation

## example
    taxa <- readTaxonData("fossil_ranges.tsv")

    lambda   ~ dnExponential(10)
    mu       ~ dnExponential(10)
    psi      ~ dnExponential(10)
    lambda_a ~ dnExponential(10)
    beta     ~ dnBeta(1, 1)

    origin ~ dnUniform(60, 80)
    moves.append( mvSlide(origin, weight=2) )

    tr ~ dnFBDSP(taxa=taxa, lambda=lambda, mu=mu, psi=psi, rho=1, lambda_a=lambda_a, beta=beta, origin=origin)

    rec ~ dnFossilRecord(ranges=tr, reporting="firstlast")
    rec.clamp(taxa)

    n <- taxa.size()
    moves.append( mvFNPR(tr, weight=n) )
    moves.append( mvNodeTimeSlideUniform(tr, weight=n) )
    moves.append( mvRootTimeSlideUniform(tr, origin=origin, weight=2) )
    moves.append( mvFossilTipTimeUniform(tr, origin=origin, weight=n) )
    moves.append( mvSpeciationMode(tr, weight=n) )
    moves.append( mvStratigraphicRange(tr, weight=n) )
    moves.append( mvCollapseExpandFossilBranch(tr, origin, weight=n/5) )
    moves.append( mvCollapseExpandBifurcation(tr, weight=n/5) )
    moves.append( mvWilsonBalding(tr, weight=n) )

    sampled := tr.getSampledTree()
    monitors.append( mnFile(sampled, filename="output/sampled.trees", printgen=100) )

## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
