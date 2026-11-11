## name
mvSpeciationMode
## title
Speciation mode move
## description
This move changes which daughter continues the ancestral species at a node of a `dnFBDSP` tree. It does not change the topology or any age.
## details
The move picks an internal node and draws uniformly among its other states. These are the continuing daughters and, when `beta` > 0, symmetric speciation, where neither daughter continues. At a sampled ancestor, the other daughter began a species by anagenesis or by symmetric speciation whose sister left no sample. When both `lambda_a` and `beta` are positive, the move also switches between these possibilities.

Use this move in every `dnFBDSP` analysis. It is also available as `mvRotateNode`.
## authors
June Walker
## see_also
dnFossilizedBirthDeathSpeciation
mvCollapseExpandBifurcation
mvBuddingTopologyRedraw
## example
    tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)
    moves.append( mvSpeciationMode(tr, weight=taxa.size()) )

    # with symmetric speciation, it also moves nodes between budding and symmetric
    beta ~ dnBeta(1, 1)
    tr_sym ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, beta=beta, taxa=taxa)
    moves.append( mvSpeciationMode(tr_sym, weight=taxa.size()) )
## references
