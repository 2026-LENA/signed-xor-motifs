# Drosophila FlyWire v783 mixing-convergence check

Recorded: 2026-10-25

## Objective

This check tests whether the signed degree-preserving null result is sensitive
to the number of attempted edge swaps. The experiments used the sign-normalized
FlyWire v783 graph, `lump-inhibitory` randomization, non-induced counting, fixed
noncanonical neurotransmitter edges, and seeds 1–20 at each mixing depth.

New experiments were run at 25, 50, and 200 attempted swaps per edge. For a
fair comparison, the mixing-100 row below uses seeds 1–20 from the existing
100-null experiment. The final row retains all 100 mixing-100 seeds as the
primary null estimate.

Observed signed-XOR count: **4,547,397,320**

## Results

| Mixing attempts per edge | Seeds | Null mean | Median | Sample SD | CV | Null range | Observed / mean | Observed / maximum null |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 25 | 20 | 4,107,956.7 | 4,068,130.0 | 1,109,819.4 | 0.270 | 2,301,273–6,055,765 | 1,106.97 | 750.92 |
| 50 | 20 | 5,033,919.8 | 4,989,488.5 | 1,726,185.0 | 0.343 | 2,578,831–9,300,302 | 903.35 | 488.95 |
| 100 | 20 | 5,192,050.4 | 4,251,383.5 | 4,520,052.7 | 0.871 | 1,404,098–22,477,248 | 875.84 | 202.31 |
| 200 | 20 | 4,022,479.3 | 3,756,561.5 | 1,463,935.9 | 0.364 | 1,673,846–7,920,945 | 1,130.50 | 574.10 |
| 100 (primary full sample) | 100 | 5,056,056.3 | 4,618,116.5 | 2,636,911.0 | 0.522 | 1,404,098–22,477,248 | 899.40 | 202.31 |

The larger variance for the 20-seed mixing-100 subset is driven by seed 17,
which produced 22,477,248 motifs. It remains 202.31 times below observed.

There is no monotonic trend in either the mean or median as mixing increases.
Pearson correlations of same-seed counts between every pair of depths range
from -0.281 to 0.153, so no persistent same-seed count state is visible in
these samples. The 20-seed distributions are small and skewed, however, and
their differences should be treated as Monte Carlo uncertainty rather than as
precise estimates of equal means.

## Integrity checks

- Each new depth completed 20 nulls and retained 20 randomization logs, 20
  counting logs, and one observed-count log.
- No new log contains an error, failure, or warning.
- All 60 new null count logs report identical graph invariants: 2,782,075 rows,
  134,181 nodes, 1,610,737 excitatory connections, 1,097,314 inhibitory
  connections, and 74,024 fixed other connections.
- Swap acceptance remained high and stable: 99.02–99.03% for excitatory edges
  and 98.93–98.95% for inhibitory edges.
- The input, randomizer, and CUDA counter hashes match across every experiment.
- All temporary randomized CSV files were removed after counting.

Detailed summaries:

- [Mixing 25](drosophila_nulls_sign_mixing25/summary.md)
- [Mixing 50](drosophila_nulls_sign_mixing50/summary.md)
- [Mixing 100, primary 100-null sample](drosophila_nulls_sign/summary.md)
- [Mixing 200](drosophila_nulls_sign_mixing200/summary.md)

## Conclusion

The enormous observed-versus-null separation is robust to all tested mixing
depths: observed enrichment ranges from 875.84-fold to 1,130.50-fold in the
same-seed comparison, and no null approaches observed. Mixing 100 is therefore
adequate for the qualitative enrichment conclusion under this null model.

This experiment does not by itself prove formal Markov-chain convergence or
pin down the null mean precisely. A stronger final analysis could add direct
edge-overlap/autocorrelation diagnostics or larger samples at every depth.
Those refinements are lower priority than testing a biologically structured
null that preserves cell class, neuropil, spatial locality, or community.
