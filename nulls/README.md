# nulls — degree- and sign-preserving null models

C tools that generate randomized connectomes for testing whether the
signed-XOR motif is enriched beyond chance. Each randomizer rewires the input
edge list while preserving the in/out degree of every node and the sign
(excitatory/inhibitory identity) of every edge, so motif counts on the
randomized graphs form the null distribution against which the observed count
is compared (Z-score, fold-enrichment).

| file                    | input format                          |
|-------------------------|---------------------------------------|
| `randomize_drosophila.c`| FlyWire / FAFB edge list              |
| `randomize_celegans.c`  | C. elegans (Cook et al.) connectome   |

## Build

The bundled `uthash.h` is the only dependency. From the repository root:

```bash
cc -O3 -march=native -o nulls/randomize_drosophila nulls/randomize_drosophila.c -lm
cc -O3 -march=native -o nulls/randomize_celegans nulls/randomize_celegans.c -lm
```

## Drosophila workflow

### 1. Normalize the graph

For a strictly binary signed graph, collapse FlyWire rows by source, target,
and sign. GABA and GLUT become one canonical inhibitory class:

```bash
python3 nulls/normalize_drosophila.py \
  data/drosophila/flywire_fafb_v783_connections.csv \
  data/drosophila/flywire_fafb_v783_connections_sign_normalized.csv \
  --collapse sign \
  --report motif_count/results/drosophila_sign_normalization.json
```

The normalizer sums `syn_count`, writes `AGGREGATED` as the neuropil, records
input/output SHA-256 hashes, and publishes the output atomically. It uses a
temporary SQLite database so memory use stays bounded.

Use the default `--collapse nt` only when a neurotransmitter-typed multigraph
is intentional. In that representation, a directed pair carrying both GABA
and GLUT remains two inhibitory edge records to the motif counter.

### 2. Run a checkpointed null experiment

Build the CUDA counter, then run the signed-degree null:

```bash
make -C motif_count cuda

python3 nulls/run_drosophila_nulls.py \
  data/drosophila/flywire_fafb_v783_connections_sign_normalized.csv \
  --nulls 100 \
  --mixing 100 \
  --mode lump-inhibitory \
  --output-dir motif_count/results/drosophila_nulls_sign
```

`lump-inhibitory` uses one ACH layer and one combined GABA/GLUT inhibitory
layer, preserving signed in/out degree. For a neurotransmitter-typed
sensitivity analysis, pair `--collapse nt` with `--mode by-nt`; that also
preserves incoming degree separately for GABA and GLUT.

The runner:

- hashes the input, randomizer, and counter into an immutable manifest;
- counts the observed graph once;
- creates, counts, and normally deletes one randomized CSV at a time;
- checkpoints each seed in `null_counts.csv`;
- regenerates JSON and Markdown summaries after every completed seed; and
- resumes completed seeds when rerun with a larger `--nulls` value.

Do not reuse an output directory with different settings or binaries; the
manifest check will reject it.

### 3. Interpret the statistics

The empirical enrichment p-value is
`(1 + number of null counts >= observed) / (N + 1)`. Thus 10 nulls cannot
produce a value below 0.0909, 100 cannot produce one below about 0.0099, and
at least 999 are needed to resolve 0.001. Treat fold enrichment and Z-scores
from small or visibly skewed null samples as descriptive, not final
significance estimates.

The motif counter is non-induced by default. Pass `--induced` to the runner
only when extra edges among the eight matched neurons should invalidate a
match.

## Tests

```bash
python3 nulls/tests/test_drosophila_workflow.py -v
```

The tests cover normalization, sign collapsing, per-layer degree
preservation, summary parsing, checkpointing, and extension of an existing
experiment.

## Data

No connectome data is included. Supply your own edge-list CSV obtained from the
original provider under its terms.

## Third-party

`uthash.h` © Troy D. Hanson, revised BSD license (notice retained in the file
header). Not covered by this repository's CC BY-NC 4.0 license.
