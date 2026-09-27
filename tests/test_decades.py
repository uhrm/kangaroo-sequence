import math

import pytest

from kangaroo_sequence import comma_children, comma_sequence_end
from kangaroo_sequence.decades import (
    BranchChoice,
    _passage,
    block_modulus,
    branch_point,
    decade_modulus,
    decade_period,
    earliest_infinite_path,
    entry_offsets,
    entry_passage,
    immortal_entries,
    start_end,
    start_passage,
    successor_lifetimes,
    transition_table,
)

ENTRY_OFFSETS_BASE_10 = (
    0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18, 19, 20, 24, 25, 26, 27, 28, 29, 30,
    35, 36, 37, 38, 39, 40, 46, 47, 48, 49, 50, 57, 58, 59, 60, 68, 69, 70, 79, 80, 90,
)  # fmt: skip


def test_period():
    assert block_modulus(10) == 72108036000
    assert decade_period(10) == (11, 924)
    assert decade_period(3) == (4, 4)


def test_entry_offsets():
    assert entry_offsets(10) == ENTRY_OFFSETS_BASE_10
    assert entry_offsets(3) == (0, 2, 3, 6)  # §9


def test_base_3_matches_paper():
    # The table in the proof of Theorem 9.1: from 3^(4h+s) + t the sequence arrives at 3^(4h+s+1) + z or ends.
    table = {
        0: {0: 2, 2: None, 3: 0, 6: 6},
        1: {0: None, 2: 3, 3: 2, 6: 0},
        2: {0: None, 2: 6, 3: 2, 6: 3},
        3: {0: 0, 2: None, 3: 3, 6: 6},
    }
    for k in range(4, 40):
        for t, z in table[k % 4].items():
            assert entry_passage(k, t, 3).exit == z


def real_passage_matches(k: int, f: int, o: int, base: int) -> None:
    """The reduced passage through decade k agrees with the one computed with the real b^k."""
    power = base**k
    real = _passage(power, base, f, o)
    reduced = _passage(decade_modulus(k, base), base, f, o)
    assert (reduced.exit, reduced.landmine) == (real.exit, real.landmine)
    assert reduced.steps_at(power) == real.steps
    assert [(d, reduced.steps_at(power, s, sl)) for d, s, sl in reduced.branch_points] == [
        (d, s) for d, s, _ in real.branch_points
    ]


@pytest.mark.parametrize("base", [3, 4, 5, 6, 10])
def test_reduced_passages_agree_with_real(base):
    K, P = decade_period(base)
    for k in [*range(K, K + min(P, 30)), 400, 401, 402]:
        for u in entry_offsets(base):
            real_passage_matches(k, 1, u, base)
        for c in range(2, base):
            real_passage_matches(k, c, 0, base)


@pytest.mark.parametrize("base", [4, 5, 10])
def test_passages_agree_with_step_by_step(base):
    for k in range(2, 7):
        for c in range(2, base):
            n, steps, branch_points = c * base**k, 0, []
            while n < base ** (k + 1) and (children := comma_children(n, base)):
                if len(children) == 2:
                    branch_points.append((n, steps))
                n, steps = children[0], steps + 1
            passage = start_passage(k, c, base)
            assert [(branch_point(k, d, base), s) for d, s, _ in passage.branch_points] == branch_points
            assert passage.steps == steps
            if passage.exit is None:
                assert n == base ** (k + 1) - passage.landmine
            else:
                assert n == base ** (k + 1) + passage.exit


@pytest.mark.parametrize("base", [3, 4, 10])
def test_start_end_agrees_with_direct(base):
    for c in range(2, base):
        for j in [*range(2, 20), 300]:
            end = start_end(c, j, base)
            assert comma_sequence_end(c * base**j, base).last == base ** (j + end.decades) - end.landmine


def test_rows_are_bijections():
    # In each decade, the 46 entries and 8 starts map one-to-one onto the 46 exits and 8 landmines.
    landmines = tuple(range(19, 83, 9))  # 10^m - r = 99...9xy with x + y = 9 (Theorem 5.1)
    for row in transition_table(10):
        passages = [*row.entries.values(), *row.starts.values()]
        assert sorted(p.exit for p in passages if p.exit is not None) == list(ENTRY_OFFSETS_BASE_10)
        assert sorted(p.landmine for p in passages if p.landmine is not None) == list(landmines)


def test_all_successor_paths_die():
    # No cycles in the decade automaton, so every successor path from a decade >= K dies.
    assert max(successor_lifetimes(10).values()) == 44
    assert max(successor_lifetimes(3).values()) == 6  # Theorem 9.1


def test_unique_immortal_entry_per_decade():
    assert all(len(entries) == 1 for entries in immortal_entries(10).values())
    assert all(len(entries) == 1 for entries in immortal_entries(3).values())


def branch_choices(root: int, base: int, count: int, indices: bool = True) -> list[BranchChoice]:
    events = (e for e in earliest_infinite_path(root, base, indices) if isinstance(e, BranchChoice))
    return [next(events) for _ in range(count)]


def test_earliest_infinite_path_base_10_matches_paper():
    # §11 and A367620.
    choices = branch_choices(20, 10, 31)
    first = choices[0]
    assert branch_point(first.k, first.d) == 19999999918
    # Counting 20 as term 1; A367620 and the paper say a(412987860) = 19999999918, one more.
    assert first.index == 412987859
    assert "".join(str(int(c.upper)) for c in choices[:30]) == "001110011101100001100101111110"  # (11.2)
    last = choices[30]
    assert branch_point(last.k, last.d) == int("2" + "9" * 84 + "27")
    assert round(math.log10(last.index), 1) == 84.8


def test_earliest_infinite_path_base_3_matches_paper():
    # Theorem 10.1: the branch points 1, 111, 12^3 11, 12^4 11, 12^7 11, ... alternate lower and higher.
    choices = branch_choices(1, 3, 12)
    assert [c.upper for c in choices] == [False, True] * 6
    assert [c.k for c in choices] == [0, 2, 5, 6, 9, 10, 13, 14, 17, 18, 21, 22]
    assert all(c.d == 1 for c in choices)
