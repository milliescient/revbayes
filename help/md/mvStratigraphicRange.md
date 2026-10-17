## name
mvStratigraphicRange
## title
Stratigraphic range move
## description
Resamples the first and last appearances of one randomly chosen taxon in a fossilized birth-death speciation process (`dnFBDSP`).
## details
The first and last appearances are latent: the record reports the bins that contain them. A tree does not carry them, so no tree move reaches them. This move redraws both uniformly within the reported bins.

It applies to `dnFBDSP` only. A `dnFBDRP` keeps the appearances in its value, where the matrix element moves sample them. Without this move the appearances never change, and the process warns at startup.
## authors
June Walker
## see_also
dnFossilizedBirthDeathSpeciation
dnFossilRecord
mvRotateNode
## example
tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

rec ~ dnFossilRecord(ranges=tr, reporting="firstlast")
rec.clamp(taxa)

moves.append( mvStratigraphicRange(tr, weight=taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
