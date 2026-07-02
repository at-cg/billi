[![GFA Tests](https://github.com/at-cg/billi/actions/workflows/test.yml/badge.svg)](https://github.com/at-cg/billi/actions/workflows/test.yml)
[![Release](https://img.shields.io/github/v/tag/at-cg/billi?label=release)](https://github.com/at-cg/billi/tags)

# Billi

Billi is a tool for identifying bubbles in pangenome or assembly graphs, represented as [bidirected graphs](https://en.wikipedia.org/wiki/Bidirected_graph) in [GFA](https://gfa-spec.github.io/GFA-spec/GFA1.html) format. Refer to our [preprint](#citation) for details. 

Billi employs an updated definition of bubbles in bidirected graphs when compared to existing definitions, e.g., [bibubbles](https://pubmed.ncbi.nlm.nih.gov/39041615/), [ultrabubbles](https://pubmed.ncbi.nlm.nih.gov/29461862/), etc. Billi can be used to identify both cyclic and acyclic subgraphs. We give rigorous proofs in our paper and argue that the bubbles reported by Billi cannot overlap each other except when they follow a nested structure. By definition, the reported bubbles will not have tip vertices (dead ends). Billi is scalable to large pangenome graphs, including the full human pangenome graphs released by HPRC. These properties make Billi a useful alternative method for unambiguously identifying variation sites and alleles in a pangenome graph.

<p align="center">
  <img src="docs/figures/bubble_nesting.png" width="700">
  <br>
  <em> Illustration of nested panbubbles and hairpins. The three red boxes and one blue box highlight the three panbubbles and one hairpin, respectively.</em>
</p>

## Table of Contents
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Usage](#usage)
  - [Decompose](#decompose)
  - [Compact](#compact)
- [Example](#example)
- [Running the test suite](#running-the-test-suite)
- [Citation](#citation)


## **Installation**
```bash
git clone https://github.com/at-cg/billi.git
cd billi
make 
```

## **Dependencies**
- **OS:** Linux.
- **Compiler:** A C++ compiler with C++17 support
- **Build tool:** `make`
- `python3` — only needed to run the test suite, not for normal builds / usage

## Usage
### Decompose
Enumerates both panbubbles and hairpins in the input graph. Graph compaction is performed internally on the input before bubble detection:

```bash
./billi decompose -i inputgraph.gfa > out.txt 
```

### Compact
Merges long non-branching paths in the input graph into single vertices, producing a smaller, equivalent GFA:
<p align="center">
  <img src="docs/figures/compaction.png" width="500">
  <br>
  <em> Visualisation of compaction operation</em>
</p>

```bash
./billi compact -i inputgraph.gfa -o compactgraph.gfa 
```

See [docs/commands.md](docs/commands.md) for the command-line options. 

## Example
```bash
./billi decompose -i test_files/edge_cases/nested.gfa > outputfile
```
<p align="center">
  <img src="docs/figures/nested.png" width="500">
  <br>
  <em>Bandage visualisation of the nested.gfa test graph</em>
</p>

The graph (nested.gfa) contains two bubbles, one completely nested in another. 

**Expected output:**

The output format is similar to [pangene](https://github.com/lh3/pangene/tree/main)

```
CC	BB	bbID	parID	side1	side2	#alleles
CC	HP	bbID	side1	side2	#alleles
CC	AL	#hap	walk	hap_id
CC
BB	0	BB:1	<s6	<s4	-1
BB	1	-1	>s1	>s3	-1
```
- Every output starts with `CC` comment lines documenting column layout, followed by one row per panbubble/hairpin found.
- **`BB`** rows are panbubbles: `bbID` (unique ID), `parID` (`-1` if top-level, or `BB:<id>`/`HP:<id>` if nested inside another panbubble/hairpin), `side1`/`side2` (the two boundary nodes, with `<`/`>` indicating the strand each is entered on), and `#alleles`.
- **`HP`** rows are hairpins: same fields as `BB`, minus `parID`.
- **`#alleles`** is `-1` when allele walks weren't computed (e.g. the input GFA has no `W`/`P` lines). When alleles are available, each one is listed on its own `AL` row directly under the corresponding `BB`/`HP` row, and the block is terminated with a `//` line.

See the [test folder](test_files) for other test cases.

## Running the test suite
```bash
python3 src/test_gfa.py --binary ./billi --test-dir test_files --verbose
```
This is the same check run on every push/PR (see [.github/workflows/test.yml](.github/workflows/test.yml)).

## Citation
If you use Billi, please cite:

> Shreeharsha G Bhat, Daanish Mahajan, and Chirag Jain. Billi: Provably Accurate and Scalable Bubble Detection in Pangenome Graphs. *bioRxiv* (2025). https://doi.org/10.1101/2025.11.21.689636
