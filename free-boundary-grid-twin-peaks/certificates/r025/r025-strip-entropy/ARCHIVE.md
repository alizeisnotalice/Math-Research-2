# Two-row strip entropy theorem

For the top row T and first two rows R, put `M=n^2-2n` and
`b=floor(n/2)`.  The proof establishes, for all n at least 3,

```
h(T)/(n! h(R)) >= binom(M+b,b)/binom(n,b) >= (n-2)^b.
```

It follows by a gate insertion injection, classified by the independent
birth set of a completion from R.  A 3 by 3 integer enumeration checks that
the 96 constructed images are legal and distinct.  Hence even a bound of
the form `h(T) <= n! C^n h(R)` cannot hold with fixed C.  The statement does
not compare h(T) with the entire one-peak count and therefore neither proves
nor disproves the original two-peak distance problem.
