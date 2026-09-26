import pytest

from kangaroo_sequence import comma_children
from kangaroo_sequence.child_graph import ChildTree, child_graph_roots, comma_parent, explore_child_tree

# (4.6)
ROOTS_BASE_10 = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18, 19, 20, 21, 25, 26, 27, 28, 29, 30, 31, 32,
    37, 38, 39, 40, 41, 42, 43, 49, 50, 51, 52, 53, 54, 62, 63, 64, 65, 74, 75, 76, 86, 87, 98,
]  # fmt: skip


def explore_naively(root: int, base: int, max_size: int) -> ChildTree | None:
    """Visit every node of the tree one by one; ``None`` if it has more than ``max_size`` nodes."""
    size = height = branch_points = leaves = largest_leaf = 0
    stack = [(root, 0)]
    while stack:
        k, depth = stack.pop()
        size += 1
        if size > max_size:
            return None
        children = comma_children(k, base)
        if not children:
            leaves += 1
            height = max(height, depth)
            largest_leaf = max(largest_leaf, k)
        elif len(children) == 2:
            branch_points += 1
        stack.extend((child, depth + 1) for child in children)
    return ChildTree(root, size, height, branch_points, leaves, largest_leaf)


def test_roots_base_10():
    assert child_graph_roots(10) == ROOTS_BASE_10


@pytest.mark.parametrize("base", [2, 3, 5, 10, 12])
def test_parent_inverts_children(base):
    for k in range(1, 20000):
        for child in comma_children(k, base):
            assert comma_parent(child, base) == k


@pytest.mark.parametrize("base", [3, 4, 5, 6, 7, 8, 9, 10])
def test_explore_agrees_with_naive(base):
    compared = 0
    for root in child_graph_roots(base):
        tree = explore_child_tree(root, base, max_digits=8)
        if tree is not None and tree.size <= 300_000:
            assert explore_naively(root, base, max_size=300_000) == tree
            compared += 1
    assert compared > 0


def test_explore_agrees_with_naive_across_branch_points():
    # Root 21 has two branch points (at 227 and 336) and about 2 * 10^3 nodes.
    tree = explore_child_tree(21)
    assert tree.branch_points == 2
    assert explore_naively(21, 10, max_size=10**5) == tree


def test_trees_with_a_single_path_are_successor_sequences():
    # (1.2) and (1.3): the sequences starting at 3, 4, 7 never meet a branch point.
    assert explore_child_tree(3) == ChildTree(3, size=2, height=1, branch_points=0, leaves=1, largest_leaf=36)
    assert explore_child_tree(4) == ChildTree(4, 199900, 199899, 0, 1, 9999945)
    assert explore_child_tree(7) == ChildTree(7, 15, 14, 0, 1, 936)


def test_all_base_10_trees_except_20_are_finite():
    # §11: root 20 has an infinite path; the longest-lived rival, rooted at 30,
    # ends at the landmine 10^365 - 82.
    trees = {root: explore_child_tree(root, max_digits=400) for root in ROOTS_BASE_10}
    assert [root for root, tree in trees.items() if tree is None] == [20]
    assert max(tree.largest_leaf for tree in trees.values() if tree is not None) == 10**365 - 82
    assert trees[30].largest_leaf == 10**365 - 82
