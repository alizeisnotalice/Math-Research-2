# Portable finite audits

Run from this directory. The terminal-bend and rooted-word scripts construct
their outputs deterministically in this directory; compare them afterward.

```sh
python3 -B audit_terminal_bend_relocation.py > verify-terminal.stdout
python3 -B audit_rooted_word_inclusion.py > verify-inclusion.stdout
python3 -B audit_one_turn_block_transpose.py > verify-one-turn.stdout
python3 -B search_monotone_rooted_inclusions.py > verify-search.stdout
```

The first three are finite exact audits for the stated graph transformations
or cut criterion. The bounded search is explicitly a witness finder; its
negative rows are not nonexistence certificates. The block-transpose result
refutes only that specified mapping, not general one-turn or average-path
comparisons.
