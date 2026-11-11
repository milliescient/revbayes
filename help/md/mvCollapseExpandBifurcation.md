## name
mvCollapseExpandBifurcation
## title
Collapse or expand a sampled ancestor at a symmetric speciation
## description
This move switches between a sampled species ending in symmetric speciation and the same species budding twice before going extinct on a `dnFBDSP` tree.
## details
The symmetric-speciation state requires two ages to coincide, so moves that change one age at a time cannot reach it. An expansion chooses which daughter buds second. It draws that budding's age uniformly between the daughter and the first budding, then draws the extinction age uniformly between the present and the second budding. A collapse reverses this proposal.

Use this move when `beta` > 0, together with `mvSpeciationMode` and `mvCollapseExpandFossilBranch`.
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
