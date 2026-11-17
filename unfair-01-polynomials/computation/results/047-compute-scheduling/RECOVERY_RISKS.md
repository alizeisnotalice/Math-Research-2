# Local total-54 complement recovery constraints

The active producer writes a chunk gzip and per-chunk receipt before its
corresponding `events.jsonl` entry.  If a worker future fails after creating
those first two files but before the parent appends the event, do **not**
delete or overwrite the orphan.  A future recovery must independently verify
the receipt binding, compressed-file SHA, uncompressed SHA, and semantic rows
before adopting it into the event ledger; otherwise it remains unknown.

Likewise, a resumed producer currently skips keys recorded in its own event
ledger.  Before treating a resumed batch as a certified result, independently
replay and verify every recorded chunk; the producer's presence check alone is
not a semantic validation.

These are recovery/certification constraints, not claims about the active
batch's status.  The currently active run had no reported failures when this
note was written.
