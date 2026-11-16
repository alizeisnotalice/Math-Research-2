(* Exact real-domain check of the proposed n=3 full-cap envelope. *)
Clear[a, alpha, beta, target, cap];
alpha = (3 Sqrt[14] - 8)/9;
beta = 2 Sqrt[3]/3 - 4/9;
cap = Min[8/9, Max[1, (a + 4/9)^2] - 4/9,
                   Max[1, (a + 8/9)^2] - 1]/2;
pieces = {
 {0 < a <= 1/9, 0},
 {1/9 <= a <= alpha, ((a + 8/9)^2 - 1)/2},
 {alpha <= a <= 5/9, 5/18},
 {5/9 <= a <= beta, ((a + 4/9)^2 - 4/9)/2},
 {beta <= a < 1, 4/9}
};
failures = Table[
  FullSimplify[Resolve[Exists[a, p[[1]] && cap != p[[2]]], Reals]],
  {p, pieces}];
order = FullSimplify[0 < 1/9 < alpha < 5/9 < beta < 1];
boundaryChecks = FullSimplify[{
 ((alpha + 8/9)^2 - 1)/2 == 5/18,
 ((beta + 4/9)^2 - 4/9)/2 == 4/9,
 ((1/9 + 8/9)^2 - 1)/2 == 0,
 ((5/9 + 4/9)^2 - 4/9)/2 == 5/18}];
passed = TrueQ[order] && And @@ (TrueQ[# === False] & /@ failures)
          && And @@ (TrueQ /@ boundaryChecks);
out = <|"schema" -> "r075-wolfram-n3-cap-piecewise-v1",
 "version" -> $Version, "pid" -> $ProcessID,
 "domain" -> "0<a<1; b upper envelope only; no halo area computation",
 "alpha" -> ToString[alpha, InputForm], "beta" -> ToString[beta, InputForm],
 "ordered_breakpoints" -> order,
 "piece_counterexample_exists" -> (ToString[#, InputForm] & /@ failures),
 "four_boundary_identities" -> boundaryChecks, "passed" -> passed|>;
Export[Last[$ScriptCommandLine], out, "RawJSON"];
Print[ToString[out, InputForm]];
Exit[If[passed, 0, 1]];
