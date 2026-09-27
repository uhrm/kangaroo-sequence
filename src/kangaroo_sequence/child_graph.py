"""The child graph G_c of comma sequences (§4, §5 and §11 of arXiv:2401.14346).

Every node has at most one parent (Lemma 5.3) and at most two children
(Theorem 5.2), so G_c is a forest of rooted trees. The roots are the numbers
below ``base**2`` without a parent (Theorem 5.5); in base 10 there are 50 of
them (4.6). The paper reports that the tree rooted at 20 is infinite and the
other 49 are finite.
"""

import math
import sys
from dataclasses import dataclass

from kangaroo_sequence.comma import block_jump, comma_children, leading_digit


def comma_parent(n: int, base: int = 10) -> int | None:
    """Return the number of which ``n`` is a comma-child (Lemma 5.3), or ``None``."""
    f, _ = leading_digit(n, base)
    k = n - ((n - f) % base) * base - f
    if k >= 1 and n in comma_children(k, base):
        return k
    return None


def child_graph_roots(base: int = 10) -> list[int]:
    """Return the roots of the child graph, i.e. the nodes that are nobody's child."""
    return [n for n in range(1, base**2) if comma_parent(n, base) is None]


@dataclass(frozen=True)
class ChildTree:
    root: int
    size: int
    """Number of nodes in the tree."""
    height: int
    """Number of edges on the longest path from the root to a leaf."""
    branch_points: int
    """Number of nodes with two children."""
    leaves: int
    """Number of leaves, i.e. landmines reached (Theorem 5.1)."""
    largest_leaf: int


def explore_child_tree(root: int, base: int = 10, max_digits: int = 1000) -> ChildTree | None:
    """Walk the entire child tree rooted at ``root``.

    Unbranched stretches are skipped with the block jumps of equation (6.4), so
    the cost grows with the number of digits reached and the number of branch
    points, not with the number of nodes. Since distinct paths never merge,
    each node is visited at most once.

    Returns ``None`` if some node reaches ``max_digits`` digits, in which case
    the tree is not shown to be finite.
    """
    limit = base**max_digits
    size = height = branch_points = leaves = largest_leaf = 0
    stack = [(root, 0)]
    while stack:
        k, depth = stack.pop()
        if k >= limit:
            return None
        steps, k = block_jump(k, base, children=True)
        size += 1 + steps
        depth += steps
        children = comma_children(k, base)
        if not children:
            leaves += 1
            height = max(height, depth)
            largest_leaf = max(largest_leaf, k)
        elif len(children) == 2:
            branch_points += 1
        stack.extend((child, depth + 1) for child in children)
    return ChildTree(root, size, height, branch_points, leaves, largest_leaf)


def main() -> None:
    base = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    max_digits = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    roots = child_graph_roots(base)
    print(f"base {base}: {len(roots)} roots, following all paths up to {max_digits} digits")
    print(f"{'root':>6} {'leaves':>7} {'branches':>9} {'height':>12} {'largest leaf':>14}")
    for root in roots:
        tree = explore_child_tree(root, base, max_digits)
        if tree is None:
            print(f"{root:>6}  still alive at {max_digits} digits")
        else:
            print(
                f"{root:>6} {tree.leaves:>7} {tree.branch_points:>9} {_magnitude(tree.height, 10):>12}"
                f" {_below_power(tree.largest_leaf, base):>14}"
            )


def _magnitude(n: int, base: int) -> str:
    """Format ``n`` exactly if short, otherwise as ``base^e`` with a fractional exponent."""
    if n < base**6:
        return str(n)
    return f"{base}^{math.log(n, base):.2f}"


def _below_power(n: int, base: int) -> str:
    """Format ``n`` exactly if short, otherwise as ``base^m - r`` with ``base^m`` the next power above ``n``.

    Leaves are landmines, which lie less than ``base**2`` below a power of ``base``.
    """
    if n < base**6:
        return str(n)
    m, power = 1, base
    while power <= n:
        m, power = m + 1, power * base
    return f"{base}^{m} - {power - n}"


if __name__ == "__main__":
    main()
