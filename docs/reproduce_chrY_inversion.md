# Reproducing the chrY (Yq11.221) inversion example: Billi vs VG vs Pangene

This page explains how to reproduce the bubble comparison shown in the manuscript figure for the recurrent 46-kbp inversion in human chromosome Y.

**Input graph:** HPRCv1.1 chrY, available at <https://zenodo.org/records/23217792> (`chrY-v1.gfa`).

## 1. Run the tools

Set the tool locations (or put them on your `PATH`):

```bash
BILLI=billi
VG=vg
K8=k8
PANGENE=pangene.js
INPUT_GFA=chrY-v1.gfa
```

Compact the graph to avoid any non-branching paths.

```bash
$BILLI compact -i $INPUT_GFA -o chrY-v1.compact.gfa
$BILLI compact -i $INPUT_GFA -r -o chrY-v1.reverseincluded.compact.gfa
```

Billi (compute panbubbles and hairpins):

```bash
$BILLI decompose -i chrY-v1.compact.gfa -w > billi_output
```

Pangene (compute bibubbles). The option `-m 200` increases sensitivity; the default is 100.

```bash
$K8 $PANGENE call -m 200 chrY-v1.reverseincluded.compact.gfa > pangene_output
```

VG (compute snarls, and then convert to JSON).

```bash
$VG snarls -t 48 -n chrY-v1.compact.gfa > vg_output
$VG view -R vg_output > vg_output.json
```

## 2. What to look for

The hairpin starts at node 57770. To see the loop, open `chrY-v1.compact.gfa` in [Bandage](https://rrwick.github.io/Bandage/) and display the neighborhood of node 57770.

### Billi

Billi reports a hairpin at node 57770 (record `HP 2 <57770 >57770 -1`) that delineates the inverted region. It contains 50 nested panbubbles.

### VG

The snarl 57770 (forward) → 64041 (backward) is nested two levels deep:

| Level | Snarl (start → end) |
|---|---|
| Snarl | 57770 (forward) → 64041 (backward) |
| Parent | 63996 (forward) → 64009 (backward) |
| Grandparent (top-level) | 63995 (forward) → 164116 (forward) |

### Pangene

The bubble starting at node 57770 is:

```
BB  19702  21900  0  >57770  <64041  -1  2  64039,64042
```

It is nested five levels below a top-level bubble. From inner to outer, the chain is:

| Depth | Bubble ID | Boundaries |
|---|---|---|
| 6 | 19702 | `>57770` `<64041` |
| 5 | 21900 | `>64043` `>64046` |
| 4 | 21895 | `>64030` `<64037` |
| 3 | 21890 | `>64015` `<64035` |
| 2 | 21882 | `>63996` `<64009` |
| 1 (top-level) | 21881 | `>63995` `>164116` |

Each of the four innermost bubbles overlaps with a sibling bubble without nesting:

| Bubble | Overlaps with (without nesting) |
|---|---|
| `>57770` `<64041` | `>64039` `>64042` |
| `>64043` `>64046` | `<64041` `>64056` |
| `>64030` `<64037` | `>64056` `>64085` |
| `>64015` `<64035` | `>64085` `>164115` |
