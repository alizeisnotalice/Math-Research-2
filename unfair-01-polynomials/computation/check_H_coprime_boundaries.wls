h = 1 + u + u^2 + u^5 + u^6;
checks = <|"pq_sum_product_infeasible" -> (Reduce[p + q == 1 && p q == 1 && p >= 0 && q >= 0, {p,q}, Reals] === False), "a2_normalization" -> FullSimplify[Reduce[p e == 1 && p + q == 1 && 0 <= e <= 1 && p >= 0 && q >= 0, {p,q,e}, Reals] \[Equivalent] (p == 1 && q == 0 && e == 1)], "a3_p_zero" -> (Reduce[p + q == 1 && r + p q == 1 && p r == 0 && p >= 0 && q >= 0 && r >= 0 && p != 0, {p,q,r}, Reals] === False)|>;
report = <|"version" -> $Version, "checks" -> checks, "all_passed" -> And @@ Values[checks], "scope" -> "Exact real boundary systems only; full arbitrary-degree theorem is analytic."|>;
(* Exact checks print JSON; redirect output to the chosen data-volume results file. *)
Print[ExportString[report,"RawJSON"]];
Exit[If[TrueQ[report["all_passed"]],0,1]];
