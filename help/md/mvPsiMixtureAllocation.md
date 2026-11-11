## name
mvPsiMixtureAllocation
## title
Gibbs draw of a species' psi mixture component
## description
This move draws one species' fossil sampling mixture component from its full conditional distribution.
## details
When you use `psiMixtureRates`, each species in a `dnFBDRP` or `dnFBDSP` has a mixture component that scales its fossil sampling rate. The move picks a species at random and draws its component from the full conditional distribution, so every proposal is accepted. Use this move to sample these allocations, as no other move changes a species' fossil sampling mixture component.
## authors
June Walker
## see_also
dnFossilizedBirthDeathRange
dnFossilizedBirthDeathSpeciation
mvGibbsMixtureAllocation
## example
    alpha ~ dnExponential(1)
    bd ~ dnFBDRP(taxa=taxa, lambda=lambda, mu=mu, psi=psi, rho=rho, psiMixtureRates=fnDiscretizeGamma(alpha, alpha, 4))

    moves.append( mvPsiMixtureAllocation(bd, weight=taxa.size()/5) )
    components := bd.getPsiMixtureComponents()
## references
