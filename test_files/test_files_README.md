# Test Files

This folder contains GFA test cases used to validate Billi's `decompose` and `compact` subcommands. Each subfolder covers a different category of graph structure. Every subfolder includes a PDF file with the visualisations of the graphs it contains.
Each `.gfa` file is paired with a corresponding `.expected` file containing the correct output against which Billi's output is checked by the test suite (`src/test_gfa.py`).

| Folder | Description | Figures |
|--------|-------------|---------|
| [panbubble_hairpins](../test_files/panbubble_hairpins/) | Core test cases for panbubble and hairpin detection, ranging from the simplest possible bubble or hairpin up to complex nested combinations of both. | [PDF](../test_files/panbubble_hairpins/A_panbubble_hairpin_figs.pdf) |
| [allele_edgecases](../test_files/allele_edgecases/) | Test cases for allele counting and walk traversal. Graphs include `W` lines so allele counts and `AL` rows appear in Billi's output. | [PDF](../test_files/allele_edgecases/A_allele_edgecase_figs.pdf) |
| [compaction_tests](../test_files/compaction_tests/) | Test cases for the `compact` subcommand. Graphs contain unbranching linear chains that should be merged, with the expected output being the compacted GFA. | [PDF](../test_files/compaction_tests/A_compaction_tests_figs.pdf) |
| [edge_cases](../test_files/edge_cases/) | Test cases targeting structural situations designed to stress-test correctness on inputs that might cause incorrect nesting, missed bubbles, or invalid contiguity assumptions (includes example from the [paper's](https://doi.org/10.1101/2025.11.21.689636) supplementary figures). | [PDF](../test_files/edge_cases/A_edge_cases_figs.pdf) |
