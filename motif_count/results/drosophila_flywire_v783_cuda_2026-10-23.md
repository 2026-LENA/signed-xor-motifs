# Drosophila signed-XOR motif: CUDA result

Recorded: 2026-10-23

## Run

The result was produced by the CUDA counter from the repository root:

```bash
./motif_count/build/cuda/drosophila_signed_xor_motif --gpu \
  data/drosophila/flywire_fafb_v783_connections.csv
```

This was a count-only, non-induced search: neither `--out` nor `--induced`
was used. The input is the FlyWire/FAFB materialization 783 connection table.
Its repository-recorded SHA-256 is:

```text
dbe2552790e30e97746c8aa990aa4dfebc33cafb4748885cc2ee13e2319f587a  drosophila/flywire_fafb_v783_connections.csv
```

The captured console tail is preserved in
[`drosophila_flywire_v783_cuda_2026-10-23.log`](drosophila_flywire_v783_cuda_2026-10-23.log).

## Result

| Metric | Value |
|---|---:|
| Raw motif mappings | 174,409,559,947 |
| Counter-reported distinct motifs | 174,409,559,947 |
| Input rows | 3,869,878 |
| Neurons | 134,181 |
| Excitatory connection rows | 2,258,155 |
| Inhibitory connection rows | 1,519,501 |
| Other/skipped connection rows | 92,222 |
| Candidate role-6 anchors | 52,222 |
| CUDA search time | 3,061 s (51 min 1 s) |
| Wall-clock time | 51 min 2.399 s |

The last periodic progress line reports 46,067,482,624 completed GPU jobs.
That is only the last observed checkpoint, not the exact final job count and
not the motif count.

## Interpretation

The counter found **174,409,559,947 matches of the implemented signed-XOR
topological pattern** in this FlyWire graph. Each match assigns eight distinct
neurons to the motif roles and requires the twelve directed, signed edges of
the six-node XOR core plus the two-node signed-feedback extension. Feedback
variant A, variant B, or both (`AB`) can satisfy a match; an `AB` match is
counted once.

`raw` and `distinct` are equal because the excitatory `F+` and inhibitory
`F-` roles make the two arms asymmetric, so the search does not apply a
two-fold symmetry division. In this output, “distinct” is the counter's
mapping-level result; it is not evidence from a separately materialized and
deduplicated list of 174 billion neuron tuples.

## What this establishes

- The required signed-XOR wiring pattern is present in the FlyWire/FAFB v783
  graph and is extremely abundant under this counter's definition.
- The complete CUDA count finished successfully; the matching `raw` and
  `distinct` values are the expected final output format, not an error.
- The run provides a baseline count suitable for CPU/CUDA regression checks
  or comparison with null/randomized networks using exactly the same parser
  and counting options.

## Important limits

- This is a structural/topological count, not proof that these circuits
  compute XOR dynamically or are functionally active.
- Because `--induced` was not used, extra connections among the same eight
  neurons are allowed. The number is therefore a non-induced motif count.
- The parser uses source ID, target ID, and neurotransmitter type. It does not
  threshold or weight edges by `syn_count`, and it does not restrict by
  neuropil.
- Connection-row multiplicity is preserved. If the same directed neuron pair
  occurs in multiple input rows, those repeated adjacency entries can
  contribute multiplicity to the count. Thus the reported number should not
  automatically be read as 174,409,559,947 unique unordered eight-neuron
  sets.
- Without `--out`, this run does not provide motif identities, the A/B/AB
  breakdown, neuron participation frequencies, or neuropil localization.
- The absolute count is strongly affected by graph size and degree. Its main
  scientific value comes from normalized comparisons against appropriate
  degree-, sign-, and identity-preserving null models, not from the absolute
  number alone.
