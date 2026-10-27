## name
mvPsiMixtureAllocation
## title
Gibbs sampling of each taxon's psi mixture component
## description
Draws one taxon's fossil sampling mixture component from its full conditional distribution.
## details
Under `dnFBDRP` or `dnFBDSP` with `psiMixtureRates`, each species carries a latent mixture component that scales its fossil sampling rate. This is the only move that changes it. Each proposal picks a taxon at random and draws its component from its full conditional, so every proposal is accepted.
## authors
June Walker
## see_also
dnFossilizedBirthDeathRange
dnFossilizedBirthDeathSpeciation
mvGibbsMixtureAllocation
## example
    alpha ~ dnExponential(1)
    bd ~ dnFBDRP(taxa=taxa, lambda=lambda, mu=mu, psi=psi, rho=rho,
                 psiMixtureRates=fnDiscretizeGamma(alpha, alpha, 4))

    moves.append( mvPsiMixtureAllocation(bd, weight=taxa.size()) )
    components := bd.getPsiMixtureComponents()
## references
