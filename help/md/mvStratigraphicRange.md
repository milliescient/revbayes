## name
mvStratigraphicRange
## title
Stratigraphic range move
## description
This move redraws the first and last appearances of one species in a `dnFBDSP` tree.
## details
The move picks a species and redraws both appearances uniformly within the reported occurrence ages. These appearances are not part of the tree, so no tree move changes them. Use this move to sample the appearances in a `dnFBDSP` analysis.

This move is only for `dnFBDSP`. A `dnFBDRP` keeps the first and last appearances in its range matrix, where matrix element moves sample them.
## authors
June Walker
## see_also
dnFossilizedBirthDeathSpeciation
dnFossilRecord
mvSpeciationMode
## example
    tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

    rec ~ dnFossilRecord(ranges=tr, reporting="firstlast")
    rec.clamp(taxa)

    moves.append( mvStratigraphicRange(tr, weight=taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
