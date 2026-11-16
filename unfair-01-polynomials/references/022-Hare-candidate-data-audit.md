# Hare 975-candidate public-data search, 2026-11-16

Scope: public read-only sources only.  No author contact, no reconstruction of
the candidate list, and no conclusion that a list does not exist elsewhere.

## Obtained source version

* arXiv v1 abstract page: <https://arxiv.org/abs/2307.07363>, submitted
  2023-07-14.  Its source link was fetched from
  <https://arxiv.org/e-print/2307.07363>.
* Downloaded archive: `arxiv-2307.07363-source.tar.gz`, SHA-256
  `dbc360a10155406169847e269de32a27f320f99a6cb51e495f98c582865a6824`.
  The archive has exactly `main.tex` and `Fit.eps`; it contains neither a
  candidate-list file nor code/data directory.
* Extracted `main.tex`, SHA-256
  `324e8f6c6859d9b84116cbdb3c885cb9f3d22a86c14d3f36a2cb2a5acdbb6cce`.
  Lines 28--34 identify the 975 count.  Lines 500--550 give the bounded
  degree-200/four-hour protocol and aggregate Table 5.1.  The only URL located
  in the source is the author home reference at line 758; it is not a data or
  code URL.

## Public endpoints checked

* The paper's cited author URL
  <https://uwaterloo.ca/scholar/kghare/home> and requested profile root
  <https://uwaterloo.ca/scholar/kghare> presently redirect to the University
  of Waterloo home page.  No paper-specific downloads were exposed.
* The published DOI resolves to
  <https://www.tandfonline.com/doi/full/10.1080/10586458.2025.2527786>.
  Its public bibliographic/search rendering confirms the same abstract and
  count, but direct read access received Cloudflare HTTP 403.  The standard
  public supplement landing URL
  <https://www.tandfonline.com/doi/suppl/10.1080/10586458.2025.2527786>
  also returned HTTP 403; no attachment name was guessed.
* Searches for the title with `code`, `975`, GitHub, and supplementary data
  returned article/abstract mirrors and seminar listings, not a candidate list
  or project repository.

## Result and interpretation

The 975-item list and the original code were **not obtained**.  The available
source describes a catalogue of *single candidate factors* `a(x)`: each has
endpoint coefficients one and every other coefficient classified as `0`, `1`,
or strictly between them.  It is not a list of fixed two-support pairs, and it
does not provide exact fractional values or a fixed cofactor `b(x)`.

The public paper reports 7,141,686 candidate-factor shapes through degree 15,
with aggregate remainder 975.  This report neither treats those 975 cases as
reconstructed nor imports them as an exclusion rule.
