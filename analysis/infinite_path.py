"""The infinite path in the base-10 child graph (§11, A367620).

For each decade k >= K the transition table has exactly one entry 10^k + u from
which an infinite child-graph path starts. Every infinite path passes through
each decade's entry, and each node has a single parent, so the infinite path is
unique: it is A367620, starting at the root 20. From decade K on its choices at
the branch points repeat with period P decades.

Run with: uv run python analysis/infinite_path.py [number of branch points to print]
"""

import math
import sys

from kangaroo_sequence.decades import BranchChoice, branch_point, decade_period, earliest_infinite_path, immortal_entries

BASE, ROOT = 10, 20
K, P = decade_period(BASE)
count = int(sys.argv[1]) if len(sys.argv) > 1 else 40

sizes = {len(entries) for entries in immortal_entries(BASE).values()}
print(f"Immortal entries per decade (k >= {K}): {sizes}\n")


def describe(k: int, d: int) -> str:
    if k < 2:
        return str(branch_point(k, d, BASE))
    return f"{d} 9^{k - 2} {d}{BASE - 1 - d}"


print(f"The first {count} branch points of A367620 (0 = smaller child, 1 = larger child):")
print(f"{'#':>4} {'choice':>6} {'term index':>14}  branch point")
choices = []
for event in earliest_infinite_path(ROOT, BASE):
    if isinstance(event, BranchChoice):
        choices.append(event)
        index = str(event.index) if event.index < 10**12 else f"~10^{math.log10(event.index):.2f}"
        print(f"{len(choices):>4} {int(event.upper):>6} {index:>14}  {describe(event.k, event.d)}")
        if len(choices) == count:
            break

# Periodicity: collect the choices decade by decade without computing term indices.
per_decade: dict[int, list[int]] = {}
k = None
for event in earliest_infinite_path(ROOT, BASE, indices=False):
    if isinstance(event, BranchChoice):
        per_decade.setdefault(event.k, []).append(int(event.upper))
    else:
        k = event[0]
        if k >= K + 2 * P:
            break


def pattern(first: int, last: int) -> str:
    return "".join(str(c) for k in range(first, last) for c in per_decade.get(k, []))


one_period = pattern(K, K + P)
assert one_period == pattern(K + P, K + 2 * P)
print(f"\nFrom decade {K} on, the choices repeat every {P} decades:")
print(f"  {len(one_period)} branch points per period, {one_period.count('1')} of them take the larger child.")
print(f"  choices before decade {K}: {pattern(0, K)}")
print(f"  choices in one period:")
for i in range(0, len(one_period), 100):
    print("   ", one_period[i : i + 100])
