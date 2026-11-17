# n=3 rank-biased GPU trace audit

After the earlier Fraction and host-build audit, a small Metal runtime trace
was actually executed with `n=3`, `P=17`, `R=1`, `tau=2`, and four move
steps.  Its raw terminal record, per-step trace payload, resource record,
and independent CPU replay are all archived here.

`n3-p17-r1-l4-trace.verify.json` records the independent replay result:
it checks the resample/append/move sequence and state fields at all nine
ranks, with the stated rank-selection rule.  It observes accepted and
rejected root--outside and root--nonroot proposals, as well as root--root
and nonroot proposals.  The replay is an implementation check of this small
trace only.  It is not a detailed-balance enumeration, a mixing result, or
an SMC-accuracy certificate.

The older `RUN_RECEIPT.md` describes the earlier host-only compilation stage;
this file documents the later, separately replayed runtime trace and takes
precedence on whether a GPU trace was run.
