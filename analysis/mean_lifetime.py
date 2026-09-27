"""The mean lifetime of successor sequences starting at c * b^j, compared with (b-1)(b+2) / (2(b-2)).

The mean is taken over the starts c * b^j with j in one period [K, K + P) of the
decade table. Bases with a period above the given limit are skipped.

Run with: uv run python analysis/mean_lifetime.py [largest base] [largest period]
"""

import sys
from fractions import Fraction

from kangaroo_sequence.decades import decade_period, entry_offsets, start_end

largest_base = int(sys.argv[1]) if len(sys.argv) > 1 else 10
largest_period = int(sys.argv[2]) if len(sys.argv) > 2 else 3000

print(f"{'base':>4} {'period':>7} {'entries':>7} {'mean':>8} {'formula':>8}")
for base in range(3, largest_base + 1):
    K, P = decade_period(base)
    if P > largest_period:
        print(f"{base:>4} {P:>7}  skipped")
        continue
    lifetimes = [start_end(c, j, base).decades for c in range(2, base) for j in range(K, K + P)]
    mean = Fraction(sum(lifetimes), len(lifetimes))
    formula = Fraction((base - 1) * (base + 2), 2 * (base - 2))
    assert mean == formula
    print(f"{base:>4} {P:>7} {len(entry_offsets(base)):>7} {str(mean):>8} {str(formula):>8}")
