"""Decade transition tables for comma sequences.

Decade ``k`` consists of the numbers with ``k + 1`` digits, ``[b^k, b^(k+1))``.
For ``k >= 2`` a successor path can enter decade ``k`` in only two ways: at a
term ``b^k + u`` with ``u < b^2`` (an *entry*), or by starting at ``c * b^k``
with ``2 <= c <= b - 1`` (Theorem 5.4; these are the larger children of the
branch points of Theorem 5.2). It leaves the decade at ``b^(k+1) + u'`` (an
*exit*, which is the next decade's entry) or dies at a landmine
``b^(k+1) - r`` (Theorem 5.1).

Such a *passage* through decade ``k`` depends on ``b^k`` only modulo
``L = block_modulus(b)``, the lcm of all block increments ``D`` of equation
(6.4). Within a stretch of constant leading digit the path moves in whole blocks
of size ``D`` until it gets close to the next boundary, and from there on it only
depends on the distance to that boundary. So ``b^k`` can be replaced by any
``B`` with ``B = b^k (mod L)``, ``2 b^2 <= B <= b^k``: the path takes the same
steps near every boundary and ``(b^k - B) * b / D`` fewer steps in the middle of
each stretch. As ``b^k mod L`` is eventually periodic in ``k``, so are the
passages; in base 10 the period is 924.
"""

import math
from collections.abc import Iterator
from dataclasses import dataclass
from fractions import Fraction
from functools import cache

from kangaroo_sequence.comma import block_increment, comma_children, leading_digit


@cache
def block_modulus(base: int = 10) -> int:
    """The lcm of all block increments ``D_f(x)`` of equation (6.4)."""
    return math.lcm(*(block_increment(f, x, base) for f in range(1, base) for x in range(base)))


def decade_modulus(k: int, base: int = 10) -> int:
    """The value used in place of ``b^k`` to compute passages through decade ``k``.

    This is ``b^k`` itself if it is small, and otherwise the smallest number at
    least ``2 b^2`` that is congruent to ``b^k`` modulo ``block_modulus(base)``
    (plus ``L`` so that it is never below ``L``).
    """
    if k < 2:
        raise ValueError("passages are only defined for decades k >= 2")
    modulus = block_modulus(base)
    reduced = pow(base, k, modulus) + modulus * -(-2 * base * base // modulus)
    if k > reduced.bit_length():
        return reduced
    return min(base**k, reduced)


@cache
def decade_period(base: int = 10) -> tuple[int, int]:
    """Return ``(K, P)`` such that decade ``k + P`` behaves like decade ``k`` for all ``k >= K``.

    ``decade_modulus(k + 1)`` is a function of ``decade_modulus(k)``, so the
    first repetition determines the pre-period ``K`` and the period ``P``.
    """
    first_seen: dict[int, int] = {}
    k = 2
    while (B := decade_modulus(k, base)) not in first_seen:
        first_seen[B] = k
        k += 1
    return first_seen[B], k - first_seen[B]


@dataclass(frozen=True)
class Passage:
    """The successor path from an entry or start of a decade to the next decade or a landmine.

    Step counts are exact for ``b^k = modulus``; for another decade ``k`` with the
    same ``decade_modulus`` use :meth:`steps_at`.
    """

    modulus: int
    exit: int | None
    """Offset ``u`` of the first term ``b^(k+1) + u`` of the next decade, or ``None`` if the path dies."""
    landmine: int | None
    """``r`` such that the path ends at the landmine ``b^(k+1) - r``, or ``None``."""
    steps: int
    """Number of steps from the input to the exit term or the landmine."""
    slope: Fraction
    """Additional steps per unit of ``b^k - modulus``."""
    branch_points: tuple[tuple[int, int, Fraction], ...]
    """``(d, steps, slope)`` for each branch point on the path, whose larger child is ``(d + 1) * b^k``."""

    def steps_at(self, power: int, steps: int | None = None, slope: Fraction | None = None) -> int:
        """Return the number of steps for the decade with ``b^k = power``."""
        extra = (self.slope if slope is None else slope) * (power - self.modulus)
        assert extra.denominator == 1
        return (self.steps if steps is None else steps) + extra.numerator


def _children(B: int, base: int, f: int, o: int) -> list[tuple[int | None, int]]:
    """Children of ``f * B + o`` in increasing order, as ``(leading digit, offset)``.

    A leading digit of ``None`` denotes the next decade, i.e. the child ``base * B + offset``.
    """
    step = (o % base) * base
    same = o + step + f  # the child with leading digit f, if it stays below (f + 1) * B
    if f < base - 1:
        up = same + 1  # the child with leading digit f + 1
        return [(f, same)] * (same < B) + [(f + 1, up - B)] * (up >= B)
    across = o + step + 1  # the child with leading digit 1 in the next decade
    return [(None, across - B)] * (across >= B) + [(f, same)] * (same < B)


@cache
def _passage(B: int, base: int, f: int, o: int) -> Passage:
    b2 = base * base
    steps, slope, branch_points = 0, Fraction(0), []
    new_stretch = True
    while True:
        # Jump over whole blocks, stopping b^2 below the boundary so that no branch point is skipped.
        increment = block_increment(f, o % base, base)
        m = max(0, (B - b2 - o) // increment)
        o += m * increment
        steps += m * base
        if new_stretch:
            # For the real b^k, the first jump in a stretch covers (b^k - B) / D more blocks.
            slope += Fraction(base, increment)
            new_stretch = False
        children = _children(B, base, f, o)
        if not children:
            assert f == base - 1, "landmines have leading digit b - 1 (Theorem 5.1)"
            return Passage(B, None, B - o, steps, slope, tuple(branch_points))
        if len(children) == 2:
            assert children[1] == (f + 1, 0), "the larger child of a branch point is (d + 1) * b^k (Theorem 5.2)"
            branch_points.append((f, steps, slope))
        steps += 1
        new_f, o = children[0]
        if new_f is None:
            return Passage(B, o, None, steps, slope, tuple(branch_points))
        new_stretch = new_f != f
        f = new_f


def entry_passage(k: int, u: int, base: int = 10) -> Passage:
    """The passage through decade ``k`` from the term ``b^k + u``."""
    return _passage(decade_modulus(k, base), base, 1, u)


def start_passage(k: int, c: int, base: int = 10) -> Passage:
    """The passage through decade ``k`` from the term ``c * b^k``."""
    return _passage(decade_modulus(k, base), base, c, 0)


@cache
def entry_offsets(base: int = 10) -> tuple[int, ...]:
    """The offsets ``u`` for which ``b^k + u`` (``k >= 3``) is the successor of a smaller number."""
    B = 2 * base**3  # any multiple of b that is at least 2 b^2: the rule near the boundary is local
    offsets = set()
    for t in range(1, base * base + 1):
        children = _children(B, base, base - 1, B - t)
        if children and children[0][0] is None:
            offsets.add(children[0][1])
    return tuple(sorted(offsets))


@dataclass(frozen=True)
class DecadeRow:
    """All passages through one decade."""

    modulus: int
    entries: dict[int, Passage]
    """Passage from ``b^k + u`` for each ``u`` in :func:`entry_offsets`."""
    starts: dict[int, Passage]
    """Passage from ``c * b^k`` for ``2 <= c <= b - 1``."""


def decade_row(k: int, base: int = 10) -> DecadeRow:
    B = decade_modulus(k, base)
    return DecadeRow(
        B,
        {u: _passage(B, base, 1, u) for u in entry_offsets(base)},
        {c: _passage(B, base, c, 0) for c in range(2, base)},
    )


def transition_table(base: int = 10) -> list[DecadeRow]:
    """The rows for decades ``K, K + 1, ..., K + P - 1``; decade ``k >= K`` uses row ``(k - K) mod P``."""
    K, P = decade_period(base)
    return [decade_row(k, base) for k in range(K, K + P)]


def _periodic(k: int, base: int) -> int:
    """A decade at least ``K`` with the same row as decade ``k >= K``."""
    K, P = decade_period(base)
    return K + (k - K) % P


# --- Successor sequences ------------------------------------------------------------------------


@dataclass(frozen=True)
class SuccessorEnd:
    decades: int
    """``m - j`` for the landmine ``b^m - r`` reached from ``c * b^j``."""
    landmine: int
    """``r``."""


def start_end(c: int, j: int, base: int = 10) -> SuccessorEnd:
    """Where the successor sequence starting at ``c * b^j`` (``j >= 2``) ends, via the transition table."""
    passage, k = start_passage(j, c, base), j
    while passage.exit is not None:
        k += 1
        passage = entry_passage(k, passage.exit, base)
    return SuccessorEnd(k + 1 - j, passage.landmine)


def successor_lifetimes(base: int = 10) -> dict[tuple[int, int], int]:
    """For each entry state ``(k, u)`` with ``K <= k < K + P``, the number of decades until death.

    A lifetime of 1 means the path dies in decade ``k``. Raises ``ValueError``
    if some path never dies, i.e. if the successor automaton has a cycle.
    """
    K, P = decade_period(base)
    lifetimes: dict[tuple[int, int], int] = {}
    for state in ((k, u) for k in range(K, K + P) for u in entry_offsets(base)):
        path = []
        while state not in lifetimes:
            if state in path:
                raise ValueError(f"successor paths through {state} never die")
            path.append(state)
            k, u = state
            passage = entry_passage(k, u, base)
            if passage.exit is None:
                lifetimes[state] = 1
                path.pop()
                break
            state = (_periodic(k + 1, base), passage.exit)
        for earlier in reversed(path):
            k, u = earlier
            lifetimes[earlier] = 1 + lifetimes[_periodic(k + 1, base), entry_passage(k, u, base).exit]
    return lifetimes


# --- Child graph ---------------------------------------------------------------------------------


@cache
def _reach(B: int, base: int, f: int, o: int) -> frozenset[int]:
    """Exits reachable from ``f * B + o`` in the child graph, i.e. choosing freely at branch points."""
    passage = _passage(B, base, f, o)
    exits = {passage.exit} - {None}
    for d, _, _ in passage.branch_points:
        exits |= _reach(B, base, d + 1, 0)
    return frozenset(exits)


def _reach_after(passage: Passage, t: int, base: int) -> frozenset[int]:
    """Exits reachable by continuing ``passage`` past its ``t``-th branch point."""
    exits = {passage.exit} - {None}
    for d, _, _ in passage.branch_points[t + 1 :]:
        exits |= _reach(passage.modulus, base, d + 1, 0)
    return frozenset(exits)


@cache
def immortal_entries(base: int = 10) -> dict[int, frozenset[int]]:
    """For each decade ``K <= k < K + P``, the offsets ``u`` from which an infinite child-graph path starts at ``b^k + u``."""
    K, P = decade_period(base)
    alive = {k: frozenset(entry_offsets(base)) for k in range(K, K + P)}
    changed = True
    while changed:
        changed = False
        for k in range(K, K + P):
            B, after = decade_modulus(k, base), alive[_periodic(k + 1, base)]
            still = frozenset(u for u in alive[k] if _reach(B, base, 1, u) & after)
            if still != alive[k]:
                alive[k], changed = still, True
    return alive


@cache
def is_immortal(k: int, u: int, base: int = 10) -> bool:
    """Whether an infinite child-graph path starts at ``b^k + u`` (``k >= 2``)."""
    K, _ = decade_period(base)
    if k >= K:
        return u in immortal_entries(base)[_periodic(k, base)]
    return any(is_immortal(k + 1, v, base) for v in _reach(decade_modulus(k, base), base, 1, u))


def branch_point(k: int, d: int, base: int = 10) -> int:
    """The branch point in decade ``k`` with leading digit ``d`` (Theorem 5.2)."""
    if k >= 2:
        return (d + 1) * base**k - base * base + d * base + (base - 1 - d)
    return next(n for n in range(max(1, base**k), base ** (k + 1)) if leading_digit(n, base)[0] == d and len(comma_children(n, base)) == 2)


@dataclass(frozen=True)
class BranchChoice:
    k: int
    """Decade of the branch point."""
    d: int
    """Leading digit of the branch point; see :func:`branch_point`."""
    upper: bool
    """Whether the path continues with the larger child."""
    index: int | None
    """Position of the branch point in the path (the root has index 1), if requested."""


def earliest_infinite_path(root: int, base: int = 10, indices: bool = True) -> Iterator[BranchChoice | tuple[int, int]]:
    """Follow the lexicographically earliest infinite child-graph path from ``root``.

    At every branch point the smaller child is taken unless no infinite path
    continues from it. Yields a :class:`BranchChoice` for each branch point and
    ``(k, u)`` whenever the path enters decade ``k >= 2`` at ``b^k + u``.
    """
    b2 = base * base

    def reaches_immortal(n: int) -> bool:
        stack = [n]
        while stack:
            n = stack.pop()
            if n >= b2:
                if is_immortal(2, n - b2, base):
                    return True
            else:
                stack.extend(comma_children(n, base))
        return False

    # Decades 0 and 1, one term at a time.
    n, index = root, 1
    while n < b2:
        children = comma_children(n, base)
        if len(children) == 2:
            upper = not reaches_immortal(children[0])
            yield BranchChoice(0 if n < base else 1, leading_digit(n, base)[0], upper, index if indices else None)
            n = children[upper]
        elif children and reaches_immortal(children[0]):
            n = children[0]
        else:
            raise ValueError(f"no infinite path starts at {root}")
        index += 1

    k, u = 2, n - b2
    while True:
        yield k, u
        power = base**k if indices else 0
        passage, first = entry_passage(k, u, base), index
        t = 0
        while t < len(passage.branch_points):
            d, steps, slope = passage.branch_points[t]
            upper = not any(is_immortal(k + 1, v, base) for v in _reach_after(passage, t, base))
            at = first + passage.steps_at(power, steps, slope) if indices else None
            yield BranchChoice(k, d, upper, at)
            if upper:
                passage, first, t = start_passage(k, d + 1, base), at + 1 if indices else None, 0
            else:
                t += 1
        assert passage.exit is not None and is_immortal(k + 1, passage.exit, base)
        index = first + passage.steps_at(power) if indices else None
        k, u = k + 1, passage.exit
