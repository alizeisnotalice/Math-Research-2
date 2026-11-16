(* Exact scalar verification of the analytic decomposition, not a halo enumerator. *)
args = Take[$ScriptCommandLine, -2];
certificate = Import[args[[1]], "RawJSON"];
parseRational[text_String] := Module[{parts = StringSplit[text, "/"]},
  If[! AllTrue[parts, StringMatchQ[#, RegularExpression["-?[0-9]+"]] &], Exit[2]];
  If[Length[parts] == 1, FromDigits[parts[[1]]], FromDigits[parts[[1]]]/FromDigits[parts[[2]]]]
];
a = Sqrt[1/8]; u = Sqrt[999/4000]; v = Sqrt[749/4000];
s = Sqrt[499/4000]; d = 71/200; dd = 171/200; w = 999/1000; e = 1/1000;
upperPiece = 1/4 + (2 a - d) 2 a - 2 (a - d + 1/4)/2;
lowerPiece = (2 u - d) (2 u - dd) + (2 d - u - v) (2 v - dd);
crossOld = (2 u - w) (2 u - dd);
decomposition = 3/4 + e^2 + upperPiece + lowerPiece + crossOld;
expanded = 3779121/1000000 - 171 a/100 - 5273 u/1000 + 91 v/40 - 2 u v;
sw = 244823/500000 - 171 a/200 - 857 (u-v)/1000;
lo = parseRational[certificate["actual_full_plane_deletion_area"]["lower"]];
hi = parseRational[certificate["actual_full_plane_deletion_area"]["upper"]];
swlo = parseRational[certificate["corner_regions"]["SW"]["area"]["lower"]];
swhi = parseRational[certificate["corner_regions"]["SW"]["area"]["upper"]];
checks = <|
 "upper_piece_identity" -> TrueQ[FullSimplify[upperPiece == dd - 2 dd a]],
 "full_decomposition_identity" -> TrueQ[FullSimplify[decomposition == expanded]],
 "upper_lower_disjoint" -> TrueQ[dd-u > a],
 "residual_single2_below" -> TrueQ[dd-s > u],
 "cross_old_positive_width" -> TrueQ[2 u > w],
 "cross_old_separate_from_local" -> TrueQ[w-u > v],
 "cross_old_separate_from_old_single" -> TrueQ[u < w-1/4],
 "full_strict_in_producer_interval" -> TrueQ[RootReduce[expanded-lo]>0 && RootReduce[hi-expanded]>0],
 "sw_strict_in_producer_interval" -> TrueQ[RootReduce[sw-swlo]>0 && RootReduce[swhi-sw]>0],
 "sw_strict_above_half_mass" -> TrueQ[sw > 1/8],
 "full_strict_below_273over250" -> TrueQ[expanded < 273/250],
 "cross_old_corrects_initial_omission" -> TrueQ[crossOld > 0]
|>;
result = <|"status" -> If[And@@Values[checks], "PASS", "FAIL"],
 "version" -> $Version, "kernel_pid" -> $ProcessID,
 "certificate_sha256" -> IntegerString[FileHash[args[[1]], "SHA256"],16,64],
 "checks" -> checks,
 "full_formula" -> ToString[expanded, InputForm],
 "full_decimal_40_digits_display_only" -> ToString[N[expanded,40],InputForm],
 "cross_old_area_40_digits_display_only" -> ToString[N[crossOld,40],InputForm],
 "scope" -> "Exact radical identities/orderings and algebraic comparison to frozen area endpoints; no new bbox, union-area or global-cap enumeration."|>;
Export[args[[2]], result, "RawJSON"];
Print[result["status"]];
Exit[If[result["status"]=="PASS",0,1]];
