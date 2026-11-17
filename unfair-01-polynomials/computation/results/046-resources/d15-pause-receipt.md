# d15 local MPS pause receipt

At 2026-11-17 16:14 Asia/Shanghai, the only owned local d15 recovery worker
was verified as PID 93524, executing the recorded
`resume_gpu_factor_search_d15.py` command and having no child processes.

Before stopping it, an atomic durable checkpoint and its binding metadata were
copied from the independent recovery run.  The copied state is step 8000,
prefix epoch/offset 4/4, Adam counter 1, with free shape `(49152, 136)` and
theta shape `(49152, 1)`.  The matching checkpoint record is step 8000.  The
saved file SHA-256 values are:

- `state.pt`: `2ef2249755373ee2c46004e059a1af5fbffe40ff4fcbbbac40198b3cbed27586`
- `checkpoint.json`: `0eb00babe699b89c2d6df639628a4d9e945fb669d364ee4a26ec4daef226dd7c`
- `config.json`: `b634c94b88aa8285c295dba7f33fb18f788ee4f43d3ecf67f3af0dace263b7c1`
- `resume-provenance.json`: `ac187c1bcbce8bcaae52b963faa228f7744f23c202536f22ab309eaa7e22f7ab`

After the copy verification, only PID 93524 received `SIGTERM`.  Within ten
seconds it had zero RSS and VSZ and appeared as `(Python)` with state `?Es`;
its original command no longer matched the local d15 worker process search.
This is an unreaped host-process residue, not an executing MPS worker.  No
SIGKILL was sent.  Memory-pressure sampling changed from 30% to 32% free and
reported zero throttled pages.  macOS did not expose a per-process Metal-memory
counter, so zero RSS/VSZ and absence of an executing command are the recorded
local release evidence.

The checkpoint is the last known durable state.  Work submitted after its
step 8000 snapshot, if any, was deliberately not claimed as saved or complete.
No new local task was started; any continuation is reserved for the MacBook.
