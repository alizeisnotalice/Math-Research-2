# MacBook n=11 fixed-root finite audit

This package retains three completed fixed-root, one-peak SMC runs on the
11 by 11 grid, at roots 13, 24, and 25.  Each has `P=262144`, `R=128`,
`L=16`, five CPU threads, exit code zero, and zero dead particles.  The CSV,
metadata, resource, and exit records are retained.  The original start
receipts, which contain private machine identifiers, are deliberately not
archived.

`audit_remote_cpu.py` validates the exact finite ledger: seed formulas,
replicate range, finite positive normalizers, proposal totals, acceptance
ranges, and exit status.  It reports the following finite diagnostics:

| root | mean fixed-root normalizer | replicate-weight ESS |
| ---: | ---: | ---: |
| 13 | 7.8801708155e158 | 117.21 |
| 24 | 1.0682899720e160 | 115.33 |
| 25 | 5.3994841603e160 | 111.74 |

The ratios between these three roots are numerical diagnostics only.  They
do not establish the previously considered coordinate comparison, any
asymptotic root-transport estimate, or a probability statement for the
two-peak law.
