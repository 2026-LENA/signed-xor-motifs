# Local connectome datasets

The large files in this directory are local, reproducible artifacts and are
ignored by Git. They were downloaded and verified on 2026-10-07.

## Ready-to-use inputs

| Counter | Input | Rows excluding header | Status |
|---|---|---:|---|
| `signed_xor_motif_count.c` | `allen_v1/allen_v1_edges.csv` | 70,139,111 | Ready; exact five-column schema |
| `drosophila_signed_xor_motif.c` | `drosophila/flywire_fafb_v783_connections.csv` | 3,869,878 | Ready; exact FlyWire schema |
| `celegans_signed_xor_motif.c` | `celegans/cook2019_hermaphrodite_chemical_nt_approx.csv` | 3,707 | Ready, but signs are an approximation; see below |

### Allen/Billeh V1

Source: the [Allen Institute mouse V1 model](https://brain-map.org/our-research/computational-modeling/models-of-the-mouse-primary-visual-cortex),
using the recurrent V1 files from the linked
[GLIF Network Dropbox directory](https://www.dropbox.com/scl/fo/9d33eedau0iavt7afwzgn/AHb7HYRjpOM0vvt09uiVEbU/old_version/GLIF%20Network/network?rlkey=gsz56gzw37lq29mf7faxtn7b1&dl=0).

The original SONATA files are retained under `allen_v1/sonata/`. The edge HDF5
contains 70,139,111 records and opens successfully with both its first and last
data blocks readable. The flattened CSV contains 54,806,770 excitatory and
15,332,341 inhibitory rows.

To recreate the flat CSV (requires `numpy` and `h5py`):

```bash
python3 motif_count/scripts/prepare_allen_v1.py \
  data/allen_v1/sonata/v1_nodes.h5 \
  data/allen_v1/sonata/v1_node_types.csv \
  data/allen_v1/sonata/v1_v1_edges.h5 \
  data/allen_v1/sonata/v1_v1_edge_types.csv \
  data/allen_v1/allen_v1_edges.csv
```

The paper's conservative threshold can be applied at runtime:

```bash
./motif_count/build/cpu/signed_xor_motif_count \
  --min-weight 0.1 data/allen_v1/allen_v1_edges.csv
```

### Drosophila FlyWire/FAFB

Source: [FlyWire Codex downloads](https://codex.flywire.ai/api/download), FAFB
materialization `783`, proofread connections with neuropil, synapse count, and
predicted neurotransmitter. The downloaded gzip passed its integrity check and
its provider MD5 is `41206440318c77418bfc7cff1cb3e0fa`.

```bash
./motif_count/build/cpu/drosophila_signed_xor_motif \
  data/drosophila/flywire_fafb_v783_connections.csv
```

### C. elegans

Connectivity and weights come from the official
[Cook et al. (2019) SI5 workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-019-1352-7/MediaObjects/41586_2019_1352_MOESM9_ESM.xlsx).
The hermaphrodite chemical matrix contains 300 neurons, 3,707 directed
neuron-to-neuron connections, and 3,669 after the C counter removes 38
self-loops.

Important: SI5 does not include the curated per-edge `Classification` column
used for the paper's reported motif counts. The ready CSV therefore joins the
primary-neurotransmitter annotations from the public
[Fenyves et al. sign-prediction dataset](https://github.com/bank-fenyves/CeConn-SignPrediction)
and applies `ACH/GLU -> EXC`, `GABA -> INH`, everything else `OTHER`. After
self-loop removal this produces 2,908 EXC, 251 INH, and 510 OTHER connections.
The paper-exact curation instead reports 2,933 EXC, 253 INH, and 483 OTHER, so
the included CSV is suitable for running and testing the counter but not for
reproducing the published C. elegans result exactly.

To recreate the approximate input using only the Python standard library:

```bash
python3 motif_count/scripts/prepare_celegans.py \
  data/celegans/raw/Cook_et_al_2019_SI5_connectome_adjacency_matrices.xlsx \
  data/celegans/raw/Fenyves_et_al_2020_sign_predictions.csv \
  data/celegans/cook2019_hermaphrodite_chemical_nt_approx.csv
```

```bash
./motif_count/build/cpu/celegans_signed_xor_motif \
  data/celegans/cook2019_hermaphrodite_chemical_nt_approx.csv
```

## SHA-256 checksums

```text
d49dd692e59e153aa3c83f5257bfc0eff51247b86d7bb183386c6d1622c70fc9  drosophila/flywire_fafb_v783_connections.csv.gz
dbe2552790e30e97746c8aa990aa4dfebc33cafb4748885cc2ee13e2319f587a  drosophila/flywire_fafb_v783_connections.csv
86d3f24c77752b1f7df4b26cdf9583461fe7e8fc76f86d2d81b3879cdaa0e389  celegans/raw/Cook_et_al_2019_SI5_connectome_adjacency_matrices.xlsx
8748f968d759157e9e955b325d59314cc0ad8f4bf4501521f380a73a674c66ff  celegans/raw/Fenyves_et_al_2020_sign_predictions.csv
a1789e18ec81c3de8ac430dcd624e16c35ddcb7d3cac976c76120fa1121b27d0  celegans/cook2019_hermaphrodite_chemical_nt_approx.csv
7fc2cc8fb17e3282c9d864c16d294d1f8ac870f35ce1a5431aa19cc13c26f50e  allen_v1/sonata/v1_node_types.csv
84c44643db9c83a835f1bf6be6fce28ab28642891b6e95366da9cb089ed9f308  allen_v1/sonata/v1_nodes.h5
1e88bfd2ba901ca9add8cf85af3e402bc74b6ba7aaa3924935aba76518b77fdf  allen_v1/sonata/v1_v1_edge_types.csv
bc772c8416befefd6465736af160f82b5d3f8f761b2cc5e8c3edbb6a9b0176fe  allen_v1/sonata/v1_v1_edges.h5
bedd54adb4419db51ebd012dd71a8a8a2e51b3a2b1898e3bdd860148061bc10a  allen_v1/allen_v1_edges.csv
```
