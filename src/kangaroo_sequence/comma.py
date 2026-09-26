"""Comma sequences (Angelini, Branicky, Resta, Sloane, Wilson; arXiv:2401.14346).

If ``k`` and ``k'`` are consecutive terms in base ``b``, then ``k' - k`` must equal
the two-digit base-``b`` number formed from the last digit of ``k`` and the first
digit of ``k'``. The smallest such ``k'`` is chosen; if none exists the sequence ends.

The fast enumeration uses equation (6.4) of the paper: while the leading digit
``f`` stays fixed, the comma numbers are ``x_j * b + f`` with
``x_j = (x + j*f) mod b``, which is periodic in ``j`` with period ``b``. Hence
``b`` steps can be skipped at once:

    k(m*b) = k + m * b * (f + sum_{j=0}^{b-1} ((x + j*f) mod b))

Note: the paper prints this as ``k + m * (b*f + sum(...))``, which is missing a
factor ``b`` on the sum; the paper's own worked example (a block sum of 460 for
f=1, b=10) agrees with the version above.
"""

from collections.abc import Iterator
from dataclasses import dataclass


def leading_digit(k: int, base: int = 10) -> tuple[int, int]:
    """Return ``(f, p)`` where ``f`` is the leading digit of ``k`` and ``p = base**(digits-1)``."""
    p = 1
    while p * base <= k:
        p *= base
    return k // p, p


def comma_successor(k: int, base: int = 10) -> int | None:
    """Return the smallest comma-successor of ``k``, or ``None`` if there is none."""
    x = k % base
    for e in range(1, base):
        candidate = k + x * base + e
        if leading_digit(candidate, base)[0] == e:
            return candidate
    return None


def comma_sequence(start: int = 1, base: int = 10) -> Iterator[int]:
    """Yield the terms of the comma sequence one by one (the naive algorithm)."""
    k: int | None = start
    while k is not None:
        yield k
        k = comma_successor(k, base)


def block_increment(f: int, x: int, base: int = 10) -> int:
    """Increase of ``k`` over ``base`` consecutive steps with leading digit ``f``, starting at last digit ``x``."""
    return base * (f + sum((x + j * f) % base for j in range(base)))


def block_jump(k: int, base: int = 10) -> tuple[int, int]:
    """Skip as many whole blocks of ``base`` steps as possible, using equation (6.4).

    Returns ``(steps, k')``: the number of terms skipped (a multiple of ``base``,
    possibly 0) and the term reached. All skipped terms, and ``k'``, keep the
    leading digit and digit count of ``k``. Since the sequence is increasing,
    it suffices that ``k' < (f + 1) * base**(digits-1)``; this also guarantees
    that ``x*base + f`` is the only valid comma number at every skipped step.
    """
    f, p = leading_digit(k, base)
    increment = block_increment(f, k % base, base)
    m = ((f + 1) * p - 1 - k) // increment
    return m * base, k + m * increment


@dataclass(frozen=True)
class CommaSequenceEnd:
    length: int
    """Number of terms in the sequence (including the starting term)."""
    last: int
    """The final term, which has no comma-successor."""


def comma_sequence_end(start: int = 1, base: int = 10, max_digits: int | None = None) -> CommaSequenceEnd | None:
    """Find the length and final term of the comma sequence starting at ``start``.

    Alternates block jumps via (6.4) with single steps across leading-digit
    changes, so the cost grows with the number of digits of the final term
    rather than with the length of the sequence.

    If ``max_digits`` is given and the sequence reaches a term with more digits,
    ``None`` is returned (comma sequences are conjectured, but not known, to be
    finite).
    """
    k, n = start, 1
    limit = base**max_digits if max_digits is not None else None
    while True:
        if limit is not None and k >= limit:
            return None
        steps, k = block_jump(k, base)
        n += steps
        successor = comma_successor(k, base)
        if successor is None:
            return CommaSequenceEnd(length=n, last=k)
        k, n = successor, n + 1
