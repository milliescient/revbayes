## name
mvCollapseExpandBifurcation

## title
Collapse or expand a sampled ancestor at a symmetric speciation

## description
On a dnFBDSP tree with symmetric speciation, moves between a sampled species that ends at a symmetric speciation and the same species budding twice before its extinction.

## details
On a dnFBDSP tree a tip is a species' extinction. A sampled species that ends at a symmetric speciation is therefore a sampled ancestor tip at that node, with its range carried on into a symmetric node of the same age that holds the two daughters. Moves that change one age at a time can never bring those two ages together, so this move maps the state to the species budding the first daughter at that age, the second below it, and ending below that, and back.

The expansion picks which daughter buds second, and draws that budding's age uniformly between the daughter and the first budding and the extinction uniformly between the present and the second budding. The collapse is determined. Candidates are all extinct fossil tips, so their number cancels.

A bifurcation is a symmetric speciation, as in the `bifurcation` mode of `getObservedSpeciation`. The collapsed state's age moves with mvFossilTipTimeUniform or mvFossilTipTimeSlideUniform, which carry the symmetric node along with the sampled ancestor.

Use it with mvSpeciationMode and mvCollapseExpandFossilBranch.

## authors
June Walker

## see_also
mvSpeciationMode
mvCollapseExpandFossilBranch
mvFossilTipTimeUniform

## example
beta ~ dnBeta(1, 1)
tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, beta=beta, taxa=taxa)
moves.append( mvSpeciationMode(tr, weight=taxa.size()) )
moves.append( mvCollapseExpandBifurcation(tr, weight=taxa.size()/5) )

## references
