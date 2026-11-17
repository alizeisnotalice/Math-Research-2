# Native Metal queue: tasks 15--23

The frozen Swift source was SHA-256
`97918305c6c2bcc02db02c36a84bf970264a5d387c2491cf3cae4f7233a3711d` and
the remote binary was SHA-256
`8235c31db552c331e3552044347d3ba784029e7ab7f027bb2440cff23b0ea457`.
Each task used its independently persisted input hash, a 131072-mask batch,
and its predecessor's handoff manifest.  Task 19, 20, and 23 retained a
normal local complete receipt.  Transport loss after remote launch affected
tasks 15--18 and 21--22; their remote journal and final object were copied
through the same strict ControlPath and revalidated batch-by-batch by
`recover_remote_stream.py` against the local immutable manifest and input.

The final handoff is `handoff-after-task23.json`, SHA-256
`57947d9307d0eba6782b2d9aff489a40942a625b5a61d08f54c4fb3e81e55fa7`.
It records each task 0--23 at its exact terminal mask boundary.  That is a
coverage/accounting statement for this finite numerical ranking search, not a
proof of the original conjecture or a strict finite exclusion.

For tasks 15--23, the independent roots-plus-C checker retained four ranking
candidates per GPU batch.  At the production `1e-5` product tolerance the
numeric counts were:

| task | candidates | numeric negative | numeric nearunknown |
| --- | ---: | ---: | ---: |
| 15 | 2048 | 2037 | 11 |
| 16 | 2048 | 2038 | 10 |
| 17 | 2048 | 2034 | 14 |
| 18 | 2048 | 2033 | 15 |
| 19 | 1024 | 1009 | 15 |
| 20 | 1024 | 1015 | 9 |
| 21 | 2048 | 2037 | 11 |
| 22 | 2048 | 2036 | 12 |
| 23 | 1024 | 1017 | 7 |

At the looser `1e-3` diagnostic tolerance every retained candidate classified
as numeric negative.  This diagnostic does not override the `1e-5` nearunknown
records and is not a certificate.
