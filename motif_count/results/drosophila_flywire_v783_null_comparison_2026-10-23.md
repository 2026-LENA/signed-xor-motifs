# Drosophila FlyWire v783 signed-XOR null comparison

Recorded: 2026-10-23; updated: 2026-10-24

## Executive result

After removing repeated neuropil rows, the implemented signed-XOR topology is
899-fold more abundant in the observed simple signed FlyWire graph than the
mean of 100 signed-degree-preserving nulls. None of the 100 signed null counts
reached the observed count, giving a one-sided empirical enrichment p-value of
`1 / 101 = 0.00990099`.

| Representation | Rows | Observed motifs | Null mean | Null sample SD | Fold enrichment | Null range |
|---|---:|---:|---:|---:|---:|---:|
| Original connection rows | 3,869,878 | 174,409,559,947 | not run | not run | not run | not run |
| Neurotransmitter-typed | 2,830,178 | 4,853,903,564 | 3,603,128.8 | 1,514,864.4 | 1,347.14 | 1,758,298–6,204,375 |
| Simple signed graph | 2,782,075 | 4,547,397,320 | 5,056,056.3 | 2,636,911.0 | 899.40 | 1,404,098–22,477,248 |

The simple signed graph is the primary topological result. It has at most one
excitatory, one inhibitory, and one other edge for each directed neuron pair.

## Input and normalization

The source was the FlyWire/FAFB materialization 783 connection table:

```text
dbe2552790e30e97746c8aa990aa4dfebc33cafb4748885cc2ee13e2319f587a  flywire_fafb_v783_connections.csv
```

The transmitter-typed normalization collapsed rows by
`(pre_root_id, post_root_id, nt_type)`:

```text
ae734d2c7986a49a690a781b970322da5fe0ce6599e9d7c16d0e10e225ad1a0e  flywire_fafb_v783_connections_normalized.csv
```

The signed normalization collapsed rows by
`(pre_root_id, post_root_id, sign)`, mapping ACH to excitatory, GABA/GLUT to
inhibitory, and all remaining types to other:

```text
2d86d7e09c6f55b29b0da60f7f8ac9c71b2c37f1bcb85cc7ebce2bbe4e9e0db6  flywire_fafb_v783_connections_sign_normalized.csv
```

Normalization preserved the total `syn_count` of 34,153,566. The typed graph
removed 1,039,700 repeated rows. The signed graph removed 1,087,803 rows,
including 45,533 directed pairs that carried both GABA and GLUT.

The original row-level motif count was 35.93 times the typed-normalized count,
showing that neuropil row multiplicity dominated the earlier absolute number.
Collapsing transmitter types to signs reduced the typed count by a further
6.31%.

## Null models

The typed sensitivity experiment used seeds 1–10 and the primary signed
experiment used seeds 1–100. Both used 100 attempted swaps per edge,
non-induced counting, and fixed noncanonical neurotransmitter edges.

| Experiment | Normalization | Randomization layers | Preserved degree |
|---|---|---|---|
| Typed sensitivity | transmitter type | ACH, GABA, GLUT separately | in/out degree per transmitter |
| Primary signed | sign | ACH and combined GABA/GLUT | excitatory/inhibitory in/out degree |

Every completed null log reported exactly the same row, neuron, sign-edge, and
candidate-anchor totals as its observed input. All 100 signed seeds have both
randomization and counting logs, no log reports an error, and the executable
hashes still match the experiment manifest. Randomized CSVs were removed after
counting; manifests, checkpoints, summaries, and logs remain.

Executable identities:

```text
837c0e03b857c253905edbaffe08648e35d66b5445681ed79ed5cc345743bf3e  randomize_drosophila
16f08e705f42777c254a2c1737e8697bde072cf86395f59ef62615ca6e144cb2  drosophila_signed_xor_motif
```

Detailed artifacts:

- [Original row-level CUDA result](drosophila_flywire_v783_cuda_2026-10-23.md)
- [Typed normalization report](drosophila_flywire_v783_normalization_2026-10-23.json)
- [Typed null summary](drosophila_nulls/summary.md)
- [Signed normalization report](drosophila_flywire_v783_sign_normalization_2026-10-23.json)
- [Signed null summary](drosophila_nulls_sign/summary.md)

## Statistical interpretation

All 110 randomized counts were far below their corresponding observed count.
The largest typed null was 782 times smaller than typed observed; the largest
signed null was 202 times smaller than signed observed.

The typed sensitivity experiment remains a ten-null preliminary result. Its
smallest attainable one-sided enrichment p-value is

```text
p_enrichment = (1 + 0) / (10 + 1) = 0.090909...
```

The primary signed experiment reached the planned 100-null milestone. No null
equaled or exceeded observed, so the add-one estimates are

```text
p_enrichment = (1 + 0) / (100 + 1) = 0.00990099...
p_two_sided = 0.01980198...
```

The one-sided value is the smallest attainable with 100 nulls. The reported
Z-scores (3,201.8 typed and 1,722.6 signed) are descriptive and should not be
treated as well-calibrated significance statistics for these skewed null
samples. The signed distribution has a visible upper tail: four seeds exceeded
10 million motifs, and seed 17 produced the maximum of 22,477,248. That maximum
is 4.44 times the null mean but still 202.3 times below observed.

## Limits and next step

- The result is structural; it does not show that matched circuits compute XOR
  dynamically or are active in vivo.
- The searches were non-induced, so additional edges among matched neurons
  were allowed.
- The counter ignores synapse weights and neuropil after normalization.
- Degree-preserving swaps do not preserve community, cell-class, spatial, or
  neuropil organization.
- The signed empirical p-value is at the 100-null resolution floor; more nulls
  are needed to estimate a smaller tail probability precisely.

The planned 100-signed-null milestone is complete. Before a final claim,
inspect mixing convergence and test a more structured null that also preserves
biological organization such as cell class, neuropil, spatial locality, or
community. Use 999 signed nulls if resolving an empirical p-value of 0.001 is
required.
