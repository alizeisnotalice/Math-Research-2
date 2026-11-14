(* Exact symbolic checks supporting the R001 tail and strict-radius bridge. *)
ClearAll[px, py, m, r, lev, dd, kk, ss, rr];
assumptions = px > 0 && py > 0 && m > 0 && r > 0 && lev > 0;
g = m/(px py) (1 + px/(2 r)) (1 + py/(2 r));
checks = <|
  "tail_strictly_decreasing" -> FullSimplify[D[g, r] < 0, Assumptions -> assumptions],
  "tail_polynomial_equivalence" -> FullSimplify[(g < lev) == (m (2 r + px) (2 r + py) < 4 lev r^2 px py), Assumptions -> assumptions],
  "uncapped_radius_margin" -> FullSimplify[ss - (kk - 1)/dd >= 1/dd, Assumptions -> dd > 0 && kk <= dd ss],
  "capped_radius_margin" -> FullSimplify[ss - rr >= 1/dd, Assumptions -> dd > 0 && kk <= dd ss && rr <= (kk - 1)/dd],
  "reflection_single_atom_identity" -> (FullSimplify[5/2 - Sqrt[2] - (1 + (Sqrt[2] - 1)^2/2)] == 0),
  "reflection_single_atom_positive" -> FullSimplify[5/2 - Sqrt[2] > 1]
|>;
Print["Wolfram version: ", $Version];
Print["Checks: ", InputForm[checks]];
If[!And @@ Values[checks], Exit[1]];
Print["ALL_EXACT_SYMBOLIC_CHECKS_PASSED"];
Exit[0];
