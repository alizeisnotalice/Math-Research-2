(* Exact real quantifier elimination for the one-axis receiver-cell formula. *)
Clear[a,b,u,v,r,x];
distanceCheck = Resolve[
  ForAll[{a,b,u,v,r}, Implies[a <= b && u <= v && r >= 0,
    Equivalent[
      Exists[x, u <= x <= v && b-r <= x <= a+r],
      2 r >= b-a && r >= b-v && r >= u-a
    ]]], Reals];
collarCheck = Resolve[
  ForAll[{a,b,x}, Implies[a <= b && b-a < 1,
    Equivalent[b-1 <= x <= a+1,
      b-1 <= x <= b || a <= x <= a+1]]], Reals];
endpointCheck = Resolve[
  ForAll[{a,b,u,v,r}, Implies[a <= b && u <= v && r >= 0,
    Equivalent[
      ForAll[x, Implies[u <= x <= v, b-x <= r && x-a <= r]],
      b-u <= r && v-a <= r
    ]]], Reals];
out = <|"schema" -> "r071-cell-formula-wolfram-v1", "version" -> $Version,
  "process_id" -> $ProcessID,
  "minimum_enclosing_radius_interval_formula" -> distanceCheck,
  "one_axis_closed_collar_union" -> collarCheck,
  "maximum_radius_endpoints" -> endpointCheck,
  "scope" -> "Exact real quantifier elimination of three one-axis linear formulas. Two-dimensional product and halo implications require the written proof."|>;
Print[ExportString[out, "RawJSON"]];
If[TrueQ[distanceCheck && collarCheck && endpointCheck], Exit[0], Exit[1]];
