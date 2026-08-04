## name
dnFossilRecord
## title
The fossil record of a set of stratigraphic ranges
## description
The fossil record of a set of species, conditioned on a fossilized birth-death range process: the probability of the observed fossil occurrences given the species ranges that produced them.
## details
This is the observation half of the fossilized birth-death range model. The range process (`dnFBDRP` or `dnFBDSP`) describes how species ranges diversify; this distribution gives the probability of the fossil record itself.

The `complete` argument selects how sampled specimens make it into the reported record: `complete=TRUE` (every sampled occurrence is reported) or `complete=FALSE` (only the oldest and youngest occurrence of each species are reported, i.e. first/last). It also determines the support of the augmented oldest age, which belongs to the range process, so both nodes use the model given here.

The `reporting` argument replaces that choice with an estimated one: the oldest and youngest occurrence of each species are always reported, and each occurrence between them is reported independently with probability `reporting`. Leaving it unset keeps the `complete` rule, and setting it to 1 recovers `complete=TRUE`. There is no value of `reporting` that gives first/last, which reports the extremes and nothing else rather than thinning the interior at a rate. Because the extremes are always reported, the number of interior occurrences is what identifies the rate against the fossil sampling rate.

The occurrences, the fossil sampling rate, its timeline, and the augmented occurrence ages are all read from the range process rather than given as arguments here. The sampling rate also appears in the range process's non-detection term, so the two nodes must share it.

## authors
June Walker
## see_also
dnFossilizedBirthDeathRange
dnBDS
mvStratigraphicRange
## example
lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)

# the birth-death range process
bd ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

# the fossil record, clamped to the observed occurrences
rec ~ dnFossilRecord(ranges=bd, complete=false)
rec.clamp(taxa)

moves.append( mvMatrixElementScale(bd, weight=taxa.size()) )
moves.append( mvMatrixElementSlide(bd, weight=taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
