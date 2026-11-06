## name
mvRotateNode

## title
Species continuation move

## description
Moves the continuation of an ancestral species between the children of a random internal node, leaving every node age and the topology unchanged.

## details
Under budding (asymmetric) speciation exactly one child of each node continues its ancestor's species, and the others begin new species there. Which one continues is a free parameter, separate from the topology and the node ages, and this move samples it. The new continuing child is drawn uniformly among the others.

With symmetric=TRUE, for a tree with symmetric speciation (beta > 0), a node may also have no continuing child: its species ends there and both children begin new ones. The move then draws uniformly among the node's other states, one per child plus symmetric, so the proposal is its own reverse. Nodes with a sampled ancestor are left alone.

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

# with symmetric speciation, also between budding and symmetric
beta ~ dnBeta(1, 1)
tr_sym ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, beta=beta, taxa=taxa)
moves.append( mvRotateNode(tr_sym, symmetric=TRUE, weight=taxa.size()/2) )

## references
