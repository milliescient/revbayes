## name
dnBDS
## title
Birth-death process with rate shifts
## description
This distribution describes a matrix of species origination and extinction times under PyRate's birth-death with rate shifts (BDS) process. It assumes that every lineage was sampled.
## details
The matrix layout and arguments are the same as for `dnFBDRP`. Because no lineage goes unsampled, each species' birth time is its own origination time.

Set `condition="sampling"` to divide each range by its probability of leaving a fossil. With this condition and the default `coexistence=FALSE`, the density is the PyRate likelihood.

Use the same moves as for `dnFBDRP`. Add a `dnFossilRecord` node for the fossil occurrences.

This model underestimates the speciation rate when many lineages go unsampled. Use `dnFBDRP` in that case.
## authors
June Walker
## see_also
dnFossilizedBirthDeathRange
dnFossilRecord
## example
    taxa <- readTaxonData("fossil_ranges.tsv")

    lambda ~ dnExponential(1)
    mu     ~ dnExponential(1)
    psi    ~ dnExponential(1)

    bd ~ dnBDS(taxa=taxa, lambda=lambda, mu=mu, psi=psi, rho=1.0, condition="sampling")

    record ~ dnFossilRecord(ranges=bd, reporting="complete")
    record.clamp(taxa)
## references
- citation: Silvestro D, Salamin N, Antonelli A, Meyer X. 2019. Improved estimation of macroevolutionary rates from fossil data using a Bayesian framework. Paleobiology 45:546-570.
  doi: 10.1017/pab.2019.23
