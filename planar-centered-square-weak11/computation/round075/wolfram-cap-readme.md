# R75 Wolfram exact scalar check

The final runner checks the proposed n=3 capacity envelope, not any halo area. With alpha=(3sqrt(14)−8)/9 and beta=2sqrt(3)/3−4/9, it asks real quantifier elimination whether any point violates each of the five proposed min-function pieces. All five existence tests return False; the breakpoint ordering and four joining identities return True. The all-square capacity reduction itself is proved separately in the core proof and is not inferred from this symbolic computation.

Wolfram15.0.0 for Mac OS X ARM actually ran the final source, launcher PID15038 and kernel PID15058. Both were absent after exit0. Source and result hashes are in the receipt. No MATLAB was started.

There was one earlier syntax failure caused by a continued Boolean expression on a new line. That launcher returned0 despite emitting ToExpression::sntx and producing no result file. It is explicitly NOT a successful mathematical check. The initial source and original receipt are preserved; its receipt has null result hash and null kernel PID, and only launcher absence was verified. Parentheses fixed the expression, after which one complete successful run produced the final result. This distinction avoids treating process exit alone as proof.

Public replay from the topic directory:

```sh
wolframscript -file computation/round075/check_orthogonal_caps.wls /Volumes/data/CodexTemp/T2/r075-cap-piecewise-replay.json
```

The source exits according to its Boolean `passed` value; a replay should also require the output file and `passed:true`, because parsing failures can bypass that exit expression. This check covers 0<a<1; the endpoint a=1 follows by direct substitution in the displayed formula. General n, unequal outside masses and area budgets are not certified by this scalar run.
