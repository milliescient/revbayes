## name
dnBDS
## title
Birth-death process with rate shifts
## description
The birth-death-with-rateshifts process describes the distribution of a matrix of species origination and extinction times under complete lineage sampling.
## details
No lineage goes unsampled, so no range has unsampled descendants to allow for, and the kernel reduces to a speciation rate where a range begins, an extinction rate where it ends, and the three rates integrated over its duration. Nothing is dropped, so nothing collapses, and the birth time each range carries is the species' own origination. That is not true of `dnFBDRP`, where dropping the lineages that left no fossil forces an integration over them, and that integration is simultaneously its p(t) factors and the collapse of origination into divergence. `getDivergenceTimes` is therefore missing here rather than degenerate, and `getOriginationTimes` reports the column.

Nothing is conditioned for you either, so the conditioning has to be declared. `condition="time"` takes the ranges as the complete set and costs nothing. `condition="sampling"` says each range reached the data by leaving a fossil and divides by its own probability of having done so. A species known only from an extant sample takes no term either way, because it is in no fossil record. `survival` and `nTaxa` describe the process that generated the ranges rather than the record, and mean here what they mean elsewhere.

`coexistence` weights each birth by how many lineages were alive to produce it, which is the per-capita structure at the event times. It is the correct term, not an assumption, and is available here as it is in `dnFBDRP`; it defaults off because the published model omits it. The omission is harmless only where the standing diversity does not vary over the span the birth times move in, and it leaves the density unnormalized, which is why this is a partial likelihood rather than a density. Left off, with `condition="sampling"`, the defaults give the likelihood PyRate uses.

Applied to a record that cannot contain the ranges which left no fossil, a lineage leaving none is invisible rather than allowed for: its birth is never credited while the speciation rate is still charged over every duration. So where many lineages go unsampled the speciation rate is underestimated, in whichever intervals those lineages fall. Use `dnFBDRP` where that matters.
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
