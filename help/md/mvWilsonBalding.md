## name
mvWilsonBalding

## title
Wilson-Balding subtree prune-and-regraft move

## description
This move prunes and regrafts a subtree on a time tree that allows sampled-ancestor tips. It also changes the age of the attachment point.

## details
The move prunes a node with its parent, excluding the root, its children and sampled-ancestor tips. It chooses uniformly among eligible older branches and fossil tips. On a branch, it draws the attachment age uniformly. A fossil becomes a sampled ancestor. The move never regrafts above the root.

Unlike `mvFNPR`, it changes the attachment age and can create or remove a sampled ancestor. On a `dnFBDSP` tree, it also chooses the continuing daughter and, if `beta` > 0, whether speciation is symmetric. It does this on its own, with no argument. Use it with `mvFNPR` and `mvCollapseExpandFossilBranch`.

## authors
June Walker

## see_also
mvFNPR
mvCollapseExpandFossilBranch
mvSpeciationMode

## example
    tr ~ dnFBDP(originAge=origin, lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)
    moves.append( mvFNPR(tr, weight=taxa.size()) )
    moves.append( mvCollapseExpandFossilBranch(tr, origin, weight=taxa.size()/5) )
    moves.append( mvWilsonBalding(tr, weight=taxa.size()) )

## references
- citation: Wilson IJ, Balding DJ (1998). Genealogical inference from microsatellite data. Genetics 150:499-510.
- citation: Gavryushkina A, Welch D, Stadler T, Drummond AJ (2014). Bayesian inference of sampled ancestor trees for epidemiology and fossil calibration. PLoS Computational Biology 10:e1003919.
