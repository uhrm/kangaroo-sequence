import itertools

import pytest

from kangaroo_sequence import (
    CommaSequenceEnd,
    block_increment,
    block_jump,
    comma_sequence,
    comma_sequence_end,
    comma_successor,
)


def test_initial_terms():
    # (1.1)
    assert list(itertools.islice(comma_sequence(1), 9)) == [1, 12, 35, 94, 135, 186, 248, 331, 344]


def test_start_3_dies_immediately():
    assert list(comma_sequence(3)) == [3, 36]
    assert comma_successor(36) is None


def test_block_increment_matches_paper_example():
    # §6: 217 repetitions of {81, 91, 1, ..., 71}, whose sum is 460, from a(1943) = 100058.
    assert block_increment(f=1, x=8) == 460


def test_block_jump_matches_paper_example():
    # From a(1943) = 100058 the sequence stays at leading digit 1 through a(4114) = 199959.
    steps, k = block_jump(100058)
    assert (steps, k) == (2170, 100058 + 217 * 460)
    assert comma_successor(k) == 199959
    assert comma_successor(199959) == 200051


@pytest.mark.parametrize("base", [2, 3, 5, 7, 10, 12, 16])
@pytest.mark.parametrize("start", range(1, 30))
def test_block_jump_agrees_with_naive(start, base):
    terms = list(itertools.islice(comma_sequence(start, base), 3000))
    for i, k in enumerate(terms):
        steps, jumped = block_jump(k, base)
        if i + steps < len(terms):
            assert terms[i + steps] == jumped


@pytest.mark.parametrize("base", [3, 4, 5, 10])
@pytest.mark.parametrize("start", range(1, 30))
def test_end_agrees_with_naive(start, base):
    fast = comma_sequence_end(start, base, max_digits=8)
    if fast is None:
        pytest.skip("sequence too long for the naive check")
    terms = list(comma_sequence(start, base))
    assert fast == CommaSequenceEnd(length=len(terms), last=terms[-1])


@pytest.mark.parametrize(
    ("start", "length", "last"),
    [
        # (1.2) and (1.3)
        (1, 2137453, 99999945),
        (2, 194697747222394, 9999999999999918),
        (3, 2, 36),
        (4, 199900, 9999945),
        (5, 19706, 999945),
        (6, 209534289952018960, 9999999999999999936),
        (7, 15, 936),
        (8, 198104936410, 9999999999972),
    ],
)
def test_lengths_and_final_terms_from_paper(start, length, last):
    assert comma_sequence_end(start) == CommaSequenceEnd(length=length, last=last)
