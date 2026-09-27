# kangaroo-sequence

Experiments with the comma sequence from *The Comma Sequence: A Simple Sequence With
Bizarre Properties* (Angelini, Branicky, Resta, Sloane, Wilson;
[arXiv:2401.14346](https://arxiv.org/abs/2401.14346)), OEIS [A121805](https://oeis.org/A121805).

```python
>>> from kangaroo_sequence import comma_sequence_end
>>> comma_sequence_end(1)
CommaSequenceEnd(length=2137453, last=99999945)
```

Run the tests with `uv run pytest`.

## Analysis

- `uv run python -m kangaroo_sequence.child_graph`: explore the 50 child trees (§4, §11).
- `uv run python analysis/successor_ends.py`: regularities between start and end points of successor sequences.
- `uv run python analysis/infinite_path.py [n]`: the unique infinite path A367620 and its periodic branch choices.
- `uv run python analysis/export_table.py`: write the 924-periodic decade transition table to `analysis/output/`.

The outputs of these scripts are kept in `analysis/output/`.
