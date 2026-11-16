# Reproduction

```sh
python3 edge_response_n3.py > edge-response-n3.json
```

The program enumerates all \(9!\) high-to-low words.  For every nonedge it
recomputes the peak set in the augmented graph directly, independently of the
response identity, using only exact integer bit operations.
