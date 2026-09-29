"""Regularities between start and end points of successor sequences.

Above b^2 every successor sequence starts at some c * 10^j (Theorem 5.4) and ends
at a landmine 10^m - r (Theorem 5.1). Above 10^366 all of them lie in the child
tree rooted at 20. This script reports how (m - j, r) depends on (c, j).

Run with: uv run python analysis/successor_ends.py
"""

from collections import Counter
from fractions import Fraction

from kangaroo_sequence import comma_sequence_end, comma_successor
from kangaroo_sequence.decades import decade_period, decade_row, entry_offsets, start_end, successor_lifetimes

BASE = 10
LANDMINES = tuple(range(19, 83, 9))
K, P = decade_period(BASE)
STARTS = range(2, BASE)

print(f"Decade rows repeat with period P = {P} from decade K = {K} on.\n")

ends = {(c, j): start_end(c, j, BASE) for c in STARTS for j in range(2, K + 2 * P)}

# 1. Periodicity of the end point in j, and from which j on it holds.
pre_period = 1 + max(j for c in STARTS for j in range(2, K + P) if ends[c, j] != ends[c, j + P])
smallest = min(q for q in range(1, P + 1) if P % q == 0 and all(ends[c, j] == ends[c, j + q] for c in STARTS for j in range(K, K + P)))
print(f"(m - j, r) depends only on c and j mod {smallest}, for all j >= {pre_period}.")
sample = [(c, j) for c in STARTS for j in (2, 3, 7, 50, 400)]
assert all(comma_sequence_end(c * BASE**j, BASE).last == BASE ** (j + ends[c, j].decades) - ends[c, j].landmine for c, j in sample)
print("  (cross-checked against direct computation for", len(sample), "starts)\n")

# 2. Decade entries and the balance in each decade.
offsets = entry_offsets(BASE)
print(f"Sequences enter decade k >= 3 at 10^k + u for {len(offsets)} offsets u:\n  {offsets}")
def balanced(row) -> bool:
    passages = [*row.entries.values(), *row.starts.values()]
    return (sorted(p.exit for p in passages if p.exit is not None) == list(offsets)
            and sorted(p.landmine for p in passages if p.landmine is not None) == list(LANDMINES))


# Decades 3..K-1 are computed with the real 10^k, decades K..K+P-1 cover the periodic range.
balance = all(balanced(decade_row(k, BASE)) for k in range(3, K + P))
print(f"In every decade k >= 3 the {len(offsets)} entries and {len(STARTS)} starts map one-to-one onto the "
      f"{len(offsets)} exits and {len(LANDMINES)} landmines: {balance}")

# Successor sequences starting below 100 also end at landmines; find the last one to die.
small_starts = [s for s in range(1, BASE**2) if not any(comma_successor(p, BASE) == s for p in range(1, s))]
small_ends = {s: comma_sequence_end(s, BASE).last for s in small_starts}
last = max(small_starts, key=small_ends.get)
m = len(str(small_ends[last]))
print(f"The {len(small_starts)} non-successors below {BASE**2} all end below 10^{m}; the last, from {last}, "
      f"ends at 10^{m} - {BASE**m - small_ends[last]}.")
print(f"  So every landmine above 10^{m} ends exactly one sequence starting at some c * 10^j.\n")

# 3. Lifetimes.
period = [ends[c, j] for c in STARTS for j in range(K, K + P)]
lifetimes = Counter(e.decades for e in period)
print(f"Mean lifetime m - j over one period: {Fraction(sum(e.decades for e in period), len(period))}"
      f" = 1 + {len(offsets)}/{len(STARTS)}")
print(f"Longest lifetime: {max(lifetimes)} decades (from a start), "
      f"{max(successor_lifetimes(BASE).values())} decades (from a decade entry)")
print("Lifetime distribution:", sorted(lifetimes.items()), "\n")

# 4. Joint distribution of (c, r).
print("Counts of (c, r) over one period of j:")
print("     r:" + "".join(f"{r:>5}" for r in LANDMINES))
for c in STARTS:
    print(f"  c = {c}:" + "".join(f"{sum(1 for j in range(K, K + P) if ends[c, j].landmine == r):>5}" for r in LANDMINES))

# 5. The excess of 7 -> 28: 7 * 10^j dies at 10^(j+1) - 28 whenever j = 5 (mod 6).
rule = all(ends[7, j] == ends[7, 5] for j in range(5, K + 2 * P, 6))
print(f"\n7 * 10^j ends at 10^(j+1) - 28 for all j >= 5 with j = 5 (mod 6): {rule} ({ends[7, 5]})")
other = sorted({j % 84 for j in range(K, K + P) if ends[7, j] == ends[7, 5] and j % 6 != 5})
count = sum(1 for j in range(K, K + P) if ends[7, j] == ends[7, 5] and j % 6 != 5)
print(f"  It also does so for {count} other j in a period: j mod 84 in {other}")

# 6. Row and column sums of the (c, r) table, and the depletion of column 28.
counts = Counter((c, ends[c, j].landmine) for c in STARTS for j in range(K, K + P))
print("\nRow sums:   ", [sum(counts[c, r] for r in LANDMINES) for c in STARTS])
print("Column sums:", [sum(counts[c, r] for c in STARTS) for r in LANDMINES])
spread = [counts[c, r] for c in STARTS for r in LANDMINES if r != 28]
rest = [counts[c, 28] for c in STARTS if c != 7]
mean_spread = sum(spread) / len(spread)
std = (sum((n - mean_spread) ** 2 for n in spread) / len(spread)) ** 0.5
print(f"Columns r != 28: mean {mean_spread:.1f}, standard deviation {std:.1f}")
print(f"Column 28 without c = 7: {sum(rest)} ends, mean {sum(rest) / len(rest):.1f}; "
      f"c = 6 has {counts[6, 28]}, {(counts[6, 28] - sum(rest) / len(rest)) / std:+.1f} standard deviations")
