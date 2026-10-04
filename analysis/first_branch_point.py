"""The index of the first branch point 19999999918 of A367620.

Counting a(1) = 20, the paper and the OEIS entry of A367620 give
a(412987860) = 19999999918. This script counts the terms with block jumps,
checks that no earlier term has two children, and compares the length of the
whole successor sequence from 20 with A330128(20) from the OEIS.
analysis/first_branch_point.c counts the same terms one step at a time.

Run with: uv run python analysis/first_branch_point.py
"""

from kangaroo_sequence import block_jump, comma_children, comma_sequence_end

TARGET = 19999999918
A330128_20 = 19278442756937613  # length of the comma sequence starting at 20, from the OEIS b-file

n, index, earlier_branch_points = 20, 1, 0  # a(1) = 20
while n != TARGET:
    steps, n = block_jump(n, children=True)  # never skips a branch point
    index += steps
    if n == TARGET:
        break
    children = comma_children(n)
    earlier_branch_points += len(children) == 2
    n, index = children[0], index + 1
children = comma_children(TARGET)
print(f"{TARGET} is term {index} (a(1) = 20) and has {len(children)} children {children}.")
print(f"Branch points before it: {earlier_branch_points}. The next term, {children[0]}, is term {index + 1}.")

end = comma_sequence_end(20)
print(f"Length of the successor sequence from 20: {end.length}; A330128(20) = {A330128_20}; equal: {end.length == A330128_20}")
