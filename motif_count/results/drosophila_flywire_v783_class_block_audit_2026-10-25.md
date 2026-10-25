# FlyWire v783 class-block audit

| Input | SHA-256 |
|---|---|
| `/home/kovalenko/signed-xor-motifs/data/drosophila/flywire_fafb_v783_connections_sign_normalized.csv` | `2d86d7e09c6f55b29b0da60f7f8ac9c71b2c37f1bcb85cc7ebce2bbe4e9e0db6` |
| `/home/kovalenko/signed-xor-motifs/data/drosophila/flywire_fafb_v783_neuron_annotations.tsv` | `b214970b55d2fbe0853bba536fdcb9e28730f4eb7ab06f600491df795da683cd` |

The graph contains 134,181 neurons,
2,708,051 ACH/GABA/GLUT signed edges, and
74,024 fixed other-transmitter edges. The annotation
table contains 139,248 neurons.

| Field | Labels in graph | Node coverage | Signed-edge coverage | Sign/class blocks | Structurally eligible edges | Frozen edges | Largest block |
|---|---:|---:|---:|---:|---:|---:|---:|
| super_class | 10 | 99.99% | 100.00% | 162 | 100.00% | 34 | 624,771 |
| cell_class | 49 | 76.41% | 48.53% | 1,076 | 99.99% | 299 | 523,135 |
| cell_sub_class | 105 | 18.43% | 4.56% | 1,079 | 99.98% | 606 | 1,109,310 |
| cell_type | 8,837 | 99.13% | 99.56% | 441,749 | 89.97% | 271,654 | 6,887 |

"Structurally eligible" means an edge lies in a sign/source-class/target-class
block with at least two edges, two distinct sources, and two distinct targets.
It is an upper bound on edges that will actually move because duplicate-edge
and self-loop rejection can still prevent individual swaps.

Recommended class field: **`super_class`**.

The recommendation selects the first (coarsest) requested hierarchy level with
at least 90.00%
graph-node coverage and 95.00%
structurally eligible signed edges. Unannotated neurons are assigned to a
single `__UNANNOTATED__` block so signed in/out degree remains preserved for
the full graph.
