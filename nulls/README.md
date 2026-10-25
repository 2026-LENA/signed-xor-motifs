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

### 2. Audit FlyWire cell-class blocks

Download the version-locked v783 publication annotations:

```bash
curl -L --fail \
  https://raw.githubusercontent.com/flyconnectome/flywire_annotations/main/supplemental_files/Supplemental_file1_neuron_annotations.tsv \
  -o data/drosophila/flywire_fafb_v783_neuron_annotations.tsv
```

The Codex API also exposes a live `consolidated_cell_types` table, but currently
requires an authenticated Codex session. The publication export above is the
official, reproducible v783 fallback and contains the hierarchy columns used
here: `super_class`, `cell_class`, `cell_sub_class`, and `cell_type`.

Audit annotation coverage and whether sign/source-class/target-class blocks
contain enough edges and distinct endpoints for double-edge swaps:

```bash
python3 nulls/analyze_drosophila_class_blocks.py \
  data/drosophila/flywire_fafb_v783_connections_sign_normalized.csv \
  data/drosophila/flywire_fafb_v783_neuron_annotations.tsv \
  --report-json motif_count/results/drosophila_flywire_v783_class_block_audit_2026-10-25.json \
  --report-md motif_count/results/drosophila_flywire_v783_class_block_audit_2026-10-25.md
```

For the 134,181-neuron graph, `super_class` is the selected level: it annotates
134,171 neurons (99.99%), covers 99.997% of signed edges at both endpoints,
and leaves only 34 of 2,708,051 signed edges in structurally frozen blocks.
`cell_class` and `cell_sub_class` cover only 76.41% and 18.43% of graph nodes;
resolved `cell_type` creates 441,749 blocks and freezes 271,654 edges.

### 3. Run a checkpointed null experiment

Build the CUDA counter, then run the selected class-constrained signed-degree
null:

```bash
make -C motif_count cuda

python3 nulls/run_drosophila_nulls.py \
  data/drosophila/flywire_fafb_v783_connections_sign_normalized.csv \
  --nulls 100 \
  --mixing 100 \
  --mode class-block \
  --annotations data/drosophila/flywire_fafb_v783_neuron_annotations.tsv \
  --class-field super_class \
  --output-dir motif_count/results/drosophila_nulls_sign_super_class
```

`class-block` partitions edges by canonical sign, source class, and target
class. It performs `(u,v),(x,y) -> (u,y),(x,v)` swaps only inside a partition,
rejecting self-loops and duplicate signed edges. This preserves every neuron’s
signed in/out degree and, more strictly, its signed edge counts to and from
each `super_class`. Missing labels share an explicit `__UNANNOTATED__` class.
Other neurotransmitters remain fixed unless `--include-other` is passed.

For an unconstrained signed baseline, `lump-inhibitory` uses one ACH layer and
one combined GABA/GLUT inhibitory layer. For a neurotransmitter-typed
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

### 4. Interpret the statistics

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

The tests cover normalization, sign collapsing, per-layer and class-block
degree preservation, hierarchy selection, manifest hashing, summary parsing,
checkpointing, and extension of an existing experiment.

## Data

No connectome data is included. Supply your own edge-list CSV obtained from the
original provider under its terms.

## Third-party

`uthash.h` © Troy D. Hanson, revised BSD license (notice retained in the file
header). Not covered by this repository's CC BY-NC 4.0 license.
