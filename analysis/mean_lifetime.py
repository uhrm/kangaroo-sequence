"""The lifetimes of successor sequences starting at c * b^j, compared with Conjecture 17 of Dougherty-Bliss and Ter-Saakov.

The lifetime m - j of the sequence from c * b^j to the landmine b^m - r is the
number of vertices (1, u, k) of its path in their graph G'_b, i.e. the path
length |P| counted by their program. Averages are taken over the starts c * b^j
with j in one period [K, K + P) of the decade table, which correspond to the
(b - 2) L(b) paths of G'_b. Bases with a period above the given limit are skipped.

Columns: the mean of |P| and the formula (b-1)(b+2) / (2(b-2)); b/2 + 1 from the
conjecture; Var|P| / (E|P|)^2, which is 1 for an exponential distribution; and
their estimate of danger intervals survived, sum C(|P|, 2) / sum |P|.

Run with: uv run python analysis/mean_lifetime.py [largest base] [largest period]
"""

import sys
from fractions import Fraction

from kangaroo_sequence.decades import decade_period, start_end

largest_base = int(sys.argv[1]) if len(sys.argv) > 1 else 10
largest_period = int(sys.argv[2]) if len(sys.argv) > 2 else 3000


def valid_count(base: int) -> int:
    """The number of vertices (1, u, k) per value of k, as computed by valid_count in rdbliss/comma."""
    return sum(1 for u in range(1, base * base) if u <= base * ((-u) % base) + base - 1)


print(f"{'base':>4} {'period':>7} {'mean':>8} {'formula':>8} {'b/2+1':>6} {'Var/E^2':>8} {'survived':>9}")
for base in range(3, largest_base + 1):
    assert valid_count(base) == (base - 1) * (base + 2) // 2
    K, P = decade_period(base)
    if P > largest_period:
        print(f"{base:>4} {P:>7}  skipped")
        continue
    lifetimes = [start_end(c, j, base).decades for c in range(2, base) for j in range(K, K + P)]
    mean = Fraction(sum(lifetimes), len(lifetimes))
    formula = Fraction((base - 1) * (base + 2), 2 * (base - 2))
    assert mean == formula
    variance = Fraction(sum(n * n for n in lifetimes), len(lifetimes)) - mean**2
    survived = Fraction(sum(n * (n - 1) // 2 for n in lifetimes), sum(lifetimes))
    print(f"{base:>4} {P:>7} {str(mean):>8} {str(formula):>8} {base / 2 + 1:>6} "
          f"{float(variance / mean**2):>8.3f} {float(survived):>9.3f}")
