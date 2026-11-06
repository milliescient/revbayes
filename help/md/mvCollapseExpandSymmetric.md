## name
mvCollapseExpandSymmetric

## title
Collapse or expand a sampled ancestor at a symmetric speciation

## description
On an extended dnFBDSP tree with symmetric speciation, moves between a sampled species that ends at a symmetric speciation and the same species budding twice before its extinction.

## details
On an extended tree a tip is a species' extinction. A sampled species that ends at a symmetric speciation is therefore a sampled ancestor tip at that node, with its range carried on into a symmetric node of the same age that holds the two daughters. Moves that change one age at a time can never bring those two ages together, so this move maps the state to the species budding the first daughter at that age, the second below it, and ending below that, and back.

The expansion picks which daughter buds second, and draws that budding's age uniformly between the daughter and the first budding and the extinction uniformly between the present and the second budding. The collapse is determined. Candidates are all extinct fossil tips, so their number cancels.

Use it with mvSpeciationMode(symmetric=TRUE) and mvCollapseExpandFossilBranch.

## authors
June Walker

## see_also
mvSpeciationMode
mvCollapseExpandFossilBranch

## example
beta ~ dnBeta(1, 1)
tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, beta=beta, taxa=taxa)
moves.append( mvSpeciationMode(tr, symmetric=TRUE, weight=taxa.size()) )
moves.append( mvCollapseExpandSymmetric(tr, weight=taxa.size()/5) )

## references
