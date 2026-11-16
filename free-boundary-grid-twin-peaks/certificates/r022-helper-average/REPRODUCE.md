# Reproduction

```sh
python3 single_valid_slot_family.py > single-valid-slot-family.json

# The script uses exact integer graph traversal/peak tests, for n=5,...,12.
# `conditional-slot-identity.md` is a proof by the delete/insert bijection and
# needs no numerical computation.
```

The script is a deterministic direct peak check on complete free-boundary grids
for (n=5,ldots,12); it uses no sampling or floating-point decisions.
