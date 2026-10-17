## name
mvRotateNode

## title
Species continuation move

## description
Moves the continuation of an ancestral species between the children of a random internal node, leaving every node age and the topology unchanged.

## details
Under budding (asymmetric) speciation exactly one child of each node continues its ancestor's species, and the others begin new species there. Which one continues is a free parameter, separate from the topology and the node ages, and this move samples it. The new continuing child is drawn uniformly among the others.

## authors
June Walker
June Walker

## see_also
dnFossilizedBirthDeathSpeciation
mvBuddingTopologyRedraw
mvStratigraphicRange

## example
tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, timeline=timeline, taxa=taxa)

# which lineage continues the ancestral species at each speciation event
moves.append( mvRotateNode(tr, weight=taxa.size()/2) )

## references
