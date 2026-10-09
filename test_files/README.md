# Test Files

This folder contains small GFA graphs used to check that Billi produces correct output. Each subfolder targets one kind of graph structure and includes a PDF that draws every graph in it.

## Folders

| Folder | Description | Figures |
|--------|-------------|---------|
| [panbubble_hairpins](../test_files/panbubble_hairpins/) | Core test cases for panbubble and hairpin detection, ranging from the simplest possible bubble or hairpin up to complex nested combinations of both. | [PDF](../test_files/panbubble_hairpins/A_panbubble_hairpin_figs.pdf) |
| [allele_edgecases](../test_files/allele_edgecases/) | Test cases for allele counting and walk traversal. Graphs include `W` lines so allele counts and `AL` rows appear in Billi's output. | [PDF](../test_files/allele_edgecases/A_allele_edgecase_figs.pdf) |
| [compaction_tests](../test_files/compaction_tests/) | Test cases for the `compact` subcommand. Graphs contain unbranching linear chains that should be merged, with the expected output being the compacted GFA. | [PDF](../test_files/compaction_tests/A_compaction_tests_figs.pdf) |
| [edge_cases](../test_files/edge_cases/) | Test cases targeting structural situations designed to stress-test correctness on inputs that might cause incorrect nesting, missed bubbles, or invalid contiguity assumptions (includes example from the [paper's](https://doi.org/10.1101/2025.11.21.689636) supplementary figures). | [PDF](../test_files/edge_cases/A_edge_cases_figs.pdf) |

## Expected outputs

Every test graph `<name>.gfa` comes with the output Billi should produce for it. Billi has two algorithms for finding bubbles, an exact one (`-e`) and a faster default heuristic, so expected outputs are stored per algorithm:

| File | Used for |
|------|----------|
| `<name>.exact.expected` | The exact algorithm (`billi decompose -e`). Present for every graph. |
| `<name>.heuristic.expected` | The heuristic algorithm (`billi decompose`). Present **only** for graphs where the heuristic output differs from the exact output. |

When a graph has no `.heuristic.expected` file, the heuristic is expected to produce the same output as the exact algorithm, and both are checked against `<name>.exact.expected`.

An expected file that contains just `ERROR` means Billi should exit with an error for that input rather than produce output. `edge_cases/snarl_overlap` is such a case because it has no tips (see *Limitations* in the main README). 

## Running the tests

From the repository root, after building Billi:

```bash
python3 src/test_script.py --binary ./billi --test-dir test_files --verbose
```

The script runs every graph twice, once with the exact algorithm and once with the heuristic, and compares each run with the matching expected file. The comparison ignores the order of output lines. The same check runs on every push and pull request (see [.github/workflows/test.yml](../.github/workflows/test.yml)).
