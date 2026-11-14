(* Scalar cross-checks only; the geometric proof is a separate artifact. *)
ClearAll[p, e, h, u, z, w, r];
assumptions = p >= 1 && 0 < e <= 1;
w = e/(9 p);
r = p + w;
tests = <|
  "strictWitnessMarginIdentity" -> FullSimplify[
    4 p^2 + e - 4 r^2 == e/9 - 4 e^2/(81 p^2)],
  "uniformStrictWitnessMargin" -> FullSimplify[
    4 p^2 + e - 4 r^2 >= 5 e/81 > 0, assumptions],
  "nextSourceExcluded" -> FullSimplify[r < 2 p, assumptions],
  "atomBandArea" -> FullSimplify[w * 2 (p - p/2) == e/9],
  "smoothBandArea" -> FullSimplify[(w/2) * 2 (7 p/8 - 5 p/8) == e/36],
  "smoothHorizontalCapture" -> FullSimplify[
    p + 3 w/4 + h < r, assumptions && 0 < h <= w/16],
  "smoothVerticalCapture" -> FullSimplify[
    7 p/8 + h < r, assumptions && 0 < h <= w/16],
  "smoothSingleSourceExcluded" -> FullSimplify[
    5 p/8 - h > Sqrt[e]/2, assumptions && 0 < h <= e/(144 p)],
  "oldHaloShiftExcluded" -> FullSimplify[
    u - z <= 0, h > 0 && u <= -h && -h <= z <= h]
|>;
result = <|"software" -> $Version,
  "scope" -> "Exact symbolic scalar checks only; not a complete geometry formalization",
  "assumptions" -> "p>=1, 0<e<=1; smoothing h as specified per formula",
  "checks" -> tests,
  "all_passed" -> And @@ Values[tests]|>;
Print[ExportString[result, "RawJSON"]];
If[And @@ Values[tests], Exit[0], Exit[1]];
