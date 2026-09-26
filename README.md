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
