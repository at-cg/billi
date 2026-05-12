import re
import subprocess
import sys
import os
import argparse

################################################################
'''
test_gfa.py - testing pipeline script
For each .gfa file found in the test directory, this script:

1.  Runs: billi decompose -e -i <input.gfa>

2.  Checks for expected errors: 
If a corresponding .expected file exists and contains "ERROR", the test passes only 
if billi exits with a non-zero exit code. For example, test_files/edge_cases/snarl_overlap.gfa 
contains a component with zero tips, hence the expected output is an ERROR.

3.  Structural invariant checks (on successful output)
- BB lines have exactly 6 fields; HP likes have exactly 5
- bbIDs are unique across all bubbles in the output
- BB bubbles have two distinct boundary nodes (side1 != side2)
- HP hairpins have the same node for both boundaries (eg: >2, <2)
- For BB and HP: declared #alleles match actual AL line count (skip if #alleles = -1)
- Ensure each allele walk starts with one boundary and ends with the other (or vice versa)
- No two BB bubbles share the same unordered boundary pair

4. Expected output comparison
Parse both actual and expected outputs and compare. The comparison is order-insensitive, i.e.,
- The order in which bubbles are reported doesn't matter
- <s,t> and <t,s> are treated as the same bubble
- The order in which the alleles are reported doesn't matter

For reference, the output format looks like 
CC   comment/header lines               (ignored)
BB   bbID parID side1 side2 #alleles    panbubble
HP   bbID side1 side2 #alleles          hairpin
AL   #hap walk hap_id                   allele walk for preceding bubble
//                                      end of bubble block

'''

def parse_walk_nodes(walk: str):
    return re.findall(r'[<>]([^<>]+)', walk)

def strip_orientation(node: str) -> str:
    return node.lstrip('<>')

def parse_output(output: str):
    bubbles = []
    current = None
    errors = []

    for lineno, raw in enumerate(output.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith('CC'):
            continue
        if line == '//':
            if current is not None:
                bubbles.append(current)
                current = None
            continue

        parts = line.split()
        tag = parts[0]

        if tag == 'BB':
            if len(parts) != 6:
                errors.append(f"Line {lineno}: malformed BB (expected 6 fields): {line!r}")
                continue
            current = {
                'type': 'BB',
                'bbID': int(parts[1]),
                'parID': parts[2],
                'side1': strip_orientation(parts[3]),
                'side2': strip_orientation(parts[4]),
                'n_alleles': int(parts[5]),
                'alleles': [],}

        elif tag == 'HP':
            if len(parts) != 5:
                errors.append(f"Line {lineno}: malformed HP (expected 5 fields): {line!r}")
                continue
            current = {
                'type': 'HP',
                'bbID': int(parts[1]),
                'parID': None,
                'side1': strip_orientation(parts[2]),
                'side2': strip_orientation(parts[3]),
                'n_alleles': int(parts[4]),
                'alleles': [],}

        elif tag == 'FB':
            current = {
                'type': 'FB',
                'bbID': int(parts[1]),
                'parID': parts[2],
                'side1': strip_orientation(parts[3]),
                'side2': strip_orientation(parts[4]),
                'n_alleles': int(parts[5]),
                'alleles': [],}

        elif tag == 'AL':
            if current is None:
                errors.append(f"Line {lineno}: AL line outside of a bubble block: {line!r}")
                continue
            if len(parts) != 4:
                errors.append(f"Line {lineno}: malformed AL (expected 4 fields): {line!r}")
                continue
            current['alleles'].append({
                'n_hap': int(parts[1]),
                'walk': parts[2],
                'hap_id': parts[3],})

        else:
            errors.append(f"Line {lineno}: unexpected tag {parts[0]!r}: {line!r}")

    if current is not None:
        bubbles.append(current)  # Handle case where file doesn't end with //

    return bubbles, errors

def check_invariants(bubbles, errors):
    seen_ids = {}

    for b in bubbles:
        btype = b['type']
        bbID = b['bbID']

        # Check: bbIDs should be unique
        if bbID in seen_ids:
            errors.append(
                f"Duplicate bbID={bbID} (types: {seen_ids[bbID]!r} and {btype!r})" 
                )
        seen_ids[bbID] = btype

        # Check: BB/FB boundary nodes must be distinct
        if btype in ('BB', 'FB'):
            if b['side1'] == b['side2']:
                errors.append(
                    f"{btype} bbID={bbID}: side1 and side2 are the same node ({b['side1']!r})"
                )

        # Check: HP side1 and side2 should be the same node (it's a self-loop)
        if btype == 'HP':
            if b['side1'] != b['side2']:
                errors.append(
                    f"HP bbID={bbID}: side1={b['side1']!r} and side2={b['side2']!r} "
                    f"differ (expected same node with opposite orientations)"
                )

        # Check: Allele count checking for BB and HP
        if btype in ('BB', 'HP') and b['n_alleles'] != -1:
            declared = b['n_alleles']
            actual = len(b['alleles'])
            if declared != actual:
                errors.append(
                    f"{btype} bbID={bbID}: declared {declared} alleles but found {actual}"
                )

        # Check: Each walk must start with side1 or side2, and end with the other
        for al in b['alleles']:
            walk = al['walk']
            nodes = parse_walk_nodes(walk)
            if len(nodes) < 2:
                errors.append(
                    f"{btype} bbID={bbID}: walk {walk!r} has fewer than 2 nodes"
                )
                continue

            first, last = nodes[0], nodes[-1]
            valid = (
                (first == b['side1'] and last == b['side2']) or
                (first == b['side2'] and last == b['side1'])
            )
            if not valid:
                errors.append(
                    f"{btype} bbID={bbID}: walk {walk!r} starts with {first!r} "
                    f"and ends with {last!r}, expected boundaries "
                    f"({b['side1']!r}, {b['side2']!r})"
                )

    # Check: the same bubble shouldn't be reported twice
    seen_pairs = {}
    for b in bubbles:
        if b['type'] in ('BB', 'FB'):
            pair = frozenset([b['side1'], b['side2']])
            if pair in seen_pairs:
                errors.append(
                    f"Duplicate bubble boundary pair {set(pair)} "
                    f"(bbIDs: {seen_pairs[pair]} and {b['bbID']})"
                )
            seen_pairs[pair] = b['bbID']

    return errors

def bubbles_to_comparable(bubbles):
    result = set()
    for b in bubbles:
        allele_walks = frozenset(al['walk'] for al in b['alleles'])
        boundary = frozenset([b['side1'], b['side2']])
        result.add((b['type'], boundary, allele_walks))
    return result

def load_expected(gfa_path: str):
    expected_path = gfa_path.replace('.gfa', '.expected')
    if not os.path.isfile(expected_path):
        return None
    with open(expected_path, 'r') as f:
        text = f.read().strip()

    if text.startswith('ERROR'):
        message = text[len('ERROR'):].strip(': ').strip() or None
        return {'type': 'error', 'message': message}
    else:
        return {'type': 'output', 'text': text}

def run_billi(binary: str, gfa_file: str):
    return subprocess.run(
        [binary, 'decompose','-e', '-i', gfa_file],
        capture_output=True,
        text=True,)

def test_file(binary: str, gfa_file: str, verbose: bool = False):
    print(f"  Testing: {gfa_file}")
    result = run_billi(binary, gfa_file)
    expected = load_expected(gfa_file)

    if expected and expected['type'] == 'error':
        if result.returncode == 0:
            print(f"    FAIL: expected an error but billi exited successfully")
            return False
        
        if expected['message']:
            combined = (result.stderr + result.stdout).strip()
            if expected['message'].lower() not in combined.lower():
                print(f"    FAIL: expected error message {expected['message']!r}")
                print(f"      got: {combined[:200]}")
                return False
        print(f"    PASS  (correctly produced an error, exit code {result.returncode})")
        return True


    # Check: the binary didn't crash
    if result.returncode != 0:
        print(f"    FAIL: billi exited with code {result.returncode}")
        if result.stderr.strip():
            print(f"    stderr: {result.stderr.strip()}")
        return False

    if verbose:
        print("    Raw output:")
        for line in result.stdout.strip().splitlines():
            print(f"      {line}")

    # Parse output
    bubbles, errors = parse_output(result.stdout)
    errors = check_invariants(bubbles, errors)

    if errors:
        print(f"    FAIL: {len(errors)} invariant violation(s):")
        for e in errors:
            print(f"      - {e}")
        return False
    
    if expected and expected['type'] == 'output':
        expected_bubbles, expected_errors = parse_output(expected['text'])
        if expected_errors:
            print(f"    WARN: expected file has parse errors, skipping comparison")
        else:
            actual_set = bubbles_to_comparable(bubbles)
            expected_set = bubbles_to_comparable(expected_bubbles)
            if actual_set != expected_set:
                print(f"    FAIL: output differs from expected")
                missing = expected_set - actual_set
                extra = actual_set - expected_set
                if missing:
                    print(f"      Missing bubbles: {missing}")
                if extra:
                    print(f"      Extra bubbles:   {extra}")
                return False
            print(f"    Note: matches expected output")
    else:
        print(f"    Note: no .expected file, checking invariants only")

    n_bb = sum(1 for b in bubbles if b['type'] == 'BB')
    n_hp = sum(1 for b in bubbles if b['type'] == 'HP')
    n_fb = sum(1 for b in bubbles if b['type'] == 'FB')
    print(f"    PASS  (BB={n_bb}, HP={n_hp}, FB={n_fb})")
    return True

def main():
    parser = argparse.ArgumentParser(description='Test billi on GFA files')
    parser.add_argument('--binary', default='./billi', help='Path to billi binary')
    parser.add_argument('--test-dir', default='test_files', help='Root test directory')
    parser.add_argument('--verbose', action='store_true', help='Print raw output for each file')
    args = parser.parse_args()

    if not os.path.isfile(args.binary):
        print(f"ERROR: Binary not found at {args.binary}")
        sys.exit(1)

    all_passed = True
    total = 0
    passed = 0

    for root, dirs, files in os.walk(args.test_dir):
        dirs.sort()  # deterministic folder order
        gfa_files = sorted(f for f in files if f.endswith('.gfa'))
        if not gfa_files:
            continue
        print(f"\n[{root}]")
        for fname in gfa_files:
            total += 1
            ok = test_file(args.binary, os.path.join(root, fname), verbose=args.verbose)
            if ok:
                passed += 1
            else:
                all_passed = False

    print(f"\n{'='*40}")
    print(f"Results: {passed}/{total} tests passed")
    print('='*40)
    sys.exit(0 if all_passed else 1)

if __name__ == '__main__':
    main()
