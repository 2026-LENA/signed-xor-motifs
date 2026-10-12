# Signed-XOR connectomic motif

Code accompanying the preprint **"From homeostasis to credit assignment: a
signed-XOR connectomic motif for local directional error signalling"**
(M. Peña Fernández, A. González Ríos, L. Lloret Iglesias, J. Marco de Lucas;
Instituto de Física de Cantabria, IFCA, CSIC–Universidad de Cantabria, 2026) [https://doi.org/10.64898/2026.06.05.730322].

The **signed-XOR motif** is an eight-neuron, twelve-edge directed signed
circuit that extends a six-neuron XOR comparator with two feedback channels of
opposite neurotransmitter identity, converting binary mismatch into directional
error signalling (one pathway potentiation, the other depression) while
respecting Dale's principle.

This repository provides the tools used in the paper. It deliberately ships
**no connectome data** (see "Data" below).

## Components

| folder         | what it is                                                        |
|----------------|-------------------------------------------------------------------|
| `motif_count/` | C tools that enumerate the signed-XOR motif in a connectome edge list (one finder per dataset format: MICrONS/Allen, FlyWire, C. elegans). |
| `nulls/`       | C tools that generate degree- and sign-preserving randomized graphs for the enrichment null model. |
| `brian2/`      | Brian2 leaky integrate-and-fire implementation that reproduces the signed-XOR truth table and its functional validation. |

Each folder has its own README with build and usage instructions.

## Quick start

The C tools need a compiler with OpenMP; the bundled `uthash.h` is the only
dependency. For example:

```bash
cd motif_count
make cpu CC=gcc
OMP_NUM_THREADS=$(nproc) ./build/cpu/signed_xor_motif_count <your_edges.csv>
```

On Apple Silicon Macs, all three motif counters also support Metal GPU search:

```bash
cd motif_count
make metal
./build/metal/signed_xor_motif_count --gpu <your_edges.csv>
```

On Linux with an NVIDIA GPU and the CUDA Toolkit:

```bash
cd motif_count
make cuda
./build/cuda/signed_xor_motif_count --gpu <your_edges.csv>
```

On Windows, use Ubuntu on WSL 2. Install it from **PowerShell as Administrator**:

```powershell
wsl --install -d Ubuntu
```

Restart if prompted, open Ubuntu, and complete its username/password setup.
Then run these commands in the **Ubuntu terminal**, replacing the repository
and input paths with your own:

```bash
sudo apt update
sudo apt install -y build-essential python3
cd "/mnt/c/Users/YOUR_USERNAME/Desktop/signed-xor-motifs/motif_count"
make clean
make cpu CC=gcc
OMP_NUM_THREADS="$(nproc)" ./build/cpu/drosophila_signed_xor_motif \
  ../data/drosophila/flywire_fafb_v783_connections.csv
```

WSL setup follows Microsoft's
[installation guide](https://learn.microsoft.com/en-us/windows/wsl/install).
See the [Windows instructions](motif_count/README.md#windows-ubuntu-on-wsl-2)
for all three counters, path handling, regression tests, and NVIDIA CUDA
through WSL 2. Windows/WSL execution has not yet been validated for this project.

See [`motif_count/README.md`](motif_count/README.md) for GPU usage, supported
options, and CPU/GPU validation commands.

The Brian2 component is a Python package:

```bash
cd brian2
pip install -r requirements.txt
python experiments/run_all.py
```

## Data

**No connectome data is distributed in this repository.** The MICrONS,
FlyWire/FAFB, *C. elegans*, and Allen Institute V1 datasets each carry their
own terms of use; obtain them from their original providers and supply your
own edge-list CSV. The tools here operate on data the user provides.

## Citing

If you use this code, please cite the accompanying preprint:

```bibtex
@article{penafernandez2026signedxor,
  title   = {From homeostasis to credit assignment: a signed-XOR connectomic
             motif for local directional error signalling},
  author  = {Pe{\~n}a Fern{\'a}ndez, Mar{\'i}a and Gonz{\'a}lez R{\'i}os, Alejandro
             and Lloret Iglesias, Lara and Marco de Lucas, Jes{\'u}s},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {[https://doi.org/10.64898/2026.06.05.730322]}
}
```

## License

The original work in this repository is licensed under the Creative Commons
Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0); see
[`LICENSE`](LICENSE). The bundled `uthash.h` is third-party (© Troy D. Hanson,
revised BSD license) and is not covered by CC BY-NC 4.0.
