# FAFB v783 `super_class` null pilot validation

This pilot validates the class-constrained randomizer on the complete
sign-normalized FAFB v783 graph. It is an execution and invariant check, not a
null-distribution significance analysis.

| Check | Result |
|---|---:|
| Input rows | 2,782,075 |
| Signed ACH/GABA/GLUT edges | 2,708,051 |
| Fixed other-transmitter edges | 74,024 |
| Sign/source-class/target-class blocks | 162 |
| Attempted swaps | 270,805,100 |
| Accepted swaps | 261,844,932 (96.69%) |
| Blocks with no accepted swap | 16 (39 edges) |
| Observed motif count | 4,547,397,320 |
| Pilot null motif count, seed 1 | 5,892,248 |

An exact streaming comparison of the source and randomized CSVs passed:

- every neuron’s signed outgoing count to each target `super_class` matched;
- every neuron’s signed incoming count from each source `super_class` matched;
- all 74,024 other-transmitter rows remained fixed;
- the signed-edge neurotransmitter, neuropil, and synapse-count metadata
  sequence matched; and
- both files contained the same 2,782,075 rows.

The 155 MB randomized CSV was removed after verification because it is a
reproducible transient artifact. The seed, input and binary hashes, command
logs, count, and timings remain in this directory. Rerunning the checkpointed
experiment with `--keep-null-csv` reproduces it.

With only one pilot null, the empirical enrichment p-value floor is 0.5 and
must not be interpreted as statistical evidence. Extend the experiment to at
least 100 seeds for a preliminary distribution.
