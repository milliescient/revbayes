## name
dnBDS
## title
Birth-death process with rate shifts
## description
The birth-death-with-rateshifts process describes the distribution of a matrix of species origination and extinction times under complete lineage sampling.
## details
The range process under complete lineage sampling: no lineage goes unsampled, so the birth time of each range is the species' own origination, and `getOriginationTimes` reports it. `getDivergenceTimes` is not available.

`condition="time"` takes the ranges as the complete set. `condition="sampling"` divides each range by its probability of leaving a fossil. `survival` and `nTaxa` mean what they do in `dnFBDRP`.

`coexistence` defaults off, as in the published model. With `condition="sampling"`, the defaults give the PyRate likelihood.

Where many lineages go unsampled, the speciation rate is underestimated. Use `dnFBDRP` there.
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
