(* Exact finite-instantiation audit for the uniform Hd stretch proof.
   The proof below establishes the displayed identities for every d>=6;
   this script independently checks the support-membership and row formulae
   for d=6,...,80 using exact symbolic coefficients. *)
Clear[u, v, d, a, b, makeRows, tests]
makeRows[n_Integer] := Module[{s, aa, bb, tt, row},
  s = {0, 1, 2, n - 1, n};
  aa = Total[(u^# a[#, v]) & /@ s];
  bb = Total[(u^# b[#, v]) & /@ s];
  tt = Expand[aa bb];
  row[k_] := Coefficient[tt, u, k] /. {b[1, v] -> 0, b[2, v] -> 0};
  {
   MemberQ[Plus @@@ Tuples[{s, s}], n + 3],
   Expand[row[n + 1]],
   Expand[row[n + 2]],
   Expand[row[n + 3]],
   Expand[row[2]]
  }
];
tests = Table[makeRows[n], {n, 6, 80}];
Print["WolframVersion=", $Version];
Print["all_dplus3_absent=", And @@ (Not[#[[1]]] & /@ tests)];
Print["all_row_dplus1_correct=", And @@ Table[
   tests[[n - 5, 2]] === a[1, v] b[n, v] + a[2, v] b[n - 1, v],
   {n, 6, 80}]];
Print["all_row_dplus2_correct=", And @@ Table[
   tests[[n - 5, 3]] === a[2, v] b[n, v], {n, 6, 80}]];
Print["all_row_dplus3_zero=", And @@ (#[[4]] === 0 & /@ tests)];
Print["all_row2_correct=", And @@ Table[
   tests[[n - 5, 5]] === a[2, v] b[0, v], {n, 6, 80}]];
Print["d6_rows=", tests[[1]]];
