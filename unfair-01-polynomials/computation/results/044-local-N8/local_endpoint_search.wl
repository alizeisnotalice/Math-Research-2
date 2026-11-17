(* Exact semialgebraic probe for the finite endpoint lemma.  Usage: wolframscript -file ... N *)
ClearAll[run];
run[n_Integer] := Module[{aa, bb, vars, coeff, eqs, ineqs, ans, line, outdir},
  aa = Join[{1}, Array[a, n - 1], {0}];
  bb = Join[{1}, Array[b, n - 1], {0}];
  vars = Join[aa[[2 ;; n]], bb[[2 ;; n]]];
  coeff[k_] := Sum[aa[[i + 1]] bb[[k - i + 1]], {i, 0, k}];
  eqs = Table[coeff[k] == 1, {k, 1, n}];
  ineqs = Join[Thread[aa >= 0], Thread[bb >= 0], {aa[[2]] > 0, bb[[2]] > 0}];
  ans = TimeConstrained[FindInstance[And @@ Join[eqs, ineqs], vars, Reals, 1], 60, $TimedOut];
  line = ToString[<|"N" -> n, "answer" -> ans|>, InputForm];
  Print[line];
  outdir = "/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/algebraic_global/general_gateway43/local_endpoint_results";
  CreateDirectory[outdir, CreateIntermediateDirectories -> True];
  Export[FileNameJoin[{outdir, "N" <> ToString[n] <> ".txt"}], line <> "\n", "String"];
];
run[ToExpression[Last[$ScriptCommandLine]]];
