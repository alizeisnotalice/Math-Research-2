# Scope of the n=7 center CRT audit

The four completed center-prime jobs recorded here exited successfully and
their saved resource records report zero swaps.  The statement concerns those
four completed jobs only.  It does not assert that every process named
`short_tail_mod` is absent: later work may run separate noncentral jobs and
is outside this certificate.

The CRT audit itself performs read-only integer verification. Its portable
version accepts raw, Pratt, reference and source paths as explicit arguments.
Original and packaged source hashes and the I/O-only changes are recorded in
the parent IO_PROVENANCE.json. The historical pending-proof raw metadata is
preserved; its interpretation uses the proved round-29 short-gap identity.
