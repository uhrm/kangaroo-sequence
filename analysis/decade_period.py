"""Block increments, the modulus M and the period of the decade table.

Run with: uv run python analysis/decade_period.py [base]
"""

import sys

from kangaroo_sequence import block_increment
from kangaroo_sequence.decades import block_modulus, decade_period, entry_offsets

base = int(sys.argv[1]) if len(sys.argv) > 1 else 10

print(f"Block increments D_f(x) of equation (6.4) in base {base}:")
for f in range(1, base):
    print(f"  f = {f}: {sorted({block_increment(f, x, base) for x in range(base)})}")

M = block_modulus(base)
factors, n, p = {}, M, 2
while n > 1:
    while n % p == 0:
        factors[p] = factors.get(p, 0) + 1
        n //= p
    p += 1
print(f"\nM = lcm of all D_f(x) = {M} = " + " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in factors.items()))

coprime = 1
for p, e in factors.items():
    if base % p:
        coprime *= p**e


def order(m: int) -> int:
    k, x = 1, base % m
    while x != 1:
        x, k = x * base % m, k + 1
    return k


orders = {p**e: order(p**e) for p, e in factors.items() if base % p}
print(f"Orders of {base} modulo the prime powers of M coprime to {base}: {orders}")
print(f"Period of {base}^k mod M: {order(coprime) if coprime > 1 else 1}")

K, P = decade_period(base)
print(f"\nThe decade table repeats with period P = {P} from decade K = {K} on.")
print(f"Entry offsets ({len(entry_offsets(base))}): {entry_offsets(base)}")
