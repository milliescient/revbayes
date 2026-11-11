## name
mvCollapseExpandFossilBranch
## title
Collapse or expand a fossil branch
## description
This move switches a fossil between a sampled ancestor and a tip on its own branch.
## details
The move picks a fossil at random. If it is a sampled ancestor, the move expands its branch by drawing the node age uniformly between the fossil and the node above. At the root, `origin` bounds this expansion. If the fossil is a tip on its own branch, the move collapses it onto its parent node, making it a sampled ancestor.

The tree distribution must allow sampled ancestors. On a `dnFBDSP` tree, sampled ancestors arise only when `lambda_a` > 0 or `beta` > 0.
## authors
June Walker
## see_also
dnFossilizedBirthDeathSpeciation
mvWilsonBalding
mvCollapseExpandBifurcation
## example
    moves.append( mvCollapseExpandFossilBranch(tr, origin, weight=taxa.size()/5) )
## references
- citation: Gavryushkina A, Welch D, Stadler T, Drummond AJ (2014). Bayesian inference of sampled ancestor trees for epidemiology and fossil calibration. PLoS Computational Biology 10:e1003919.
