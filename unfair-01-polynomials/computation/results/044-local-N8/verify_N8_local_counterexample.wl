(* Independent exact verifier for the N=8 local endpoint counterexample.
   It reads only the saved algebraic candidate, not the search program. *)
ClearAll[a, b, coeff];
candidate = ToExpression[StringTrim[Import[
  "/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/algebraic_global/general_gateway43/local_endpoint_results/N8.txt", "String"]]];
rules = candidate["answer"][[1]];
aa = Join[{1}, Array[a, 7], {0}] /. rules;
bb = Join[{1}, Array[b, 7], {0}] /. rules;
coeff[k_] := FullSimplify[Sum[aa[[i + 1]] bb[[k - i + 1]], {i, 0, k}]];
convolutionOK = And @@ Table[TrueQ[coeff[k] == 1], {k, 0, 8}];
nonnegativeOK = And @@ (TrueQ[FullSimplify[# >= 0]] & /@ Join[aa, bb]);
firstPositiveOK = TrueQ[FullSimplify[aa[[2]] > 0]] && TrueQ[FullSimplify[bb[[2]] > 0]];
endpointsOK = aa[[1]] == 1 && bb[[1]] == 1 && aa[[9]] == 0 && bb[[9]] == 0;
Print[<|"convolution_0_to_8" -> convolutionOK,
  "all_coefficients_nonnegative" -> nonnegativeOK,
  "a1_b1_strictly_positive" -> firstPositiveOK,
  "endpoint_values" -> endpointsOK,
  "a_numeric" -> N[aa, 12], "b_numeric" -> N[bb, 12]|>];
If[!(convolutionOK && nonnegativeOK && firstPositiveOK && endpointsOK), Exit[1]];
