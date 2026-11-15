(* Independent exact prefix exclusion for A=1+s x+t x^r+x^(2r+2), r=2..8.
   This writes one InputForm record per r.  Every Reduce call is over Reals. *)
$HistoryLength = 0;
Clear[s, t, b, c, r, k, u, v];
out = FileNameJoin[{DirectoryName[$InputFileName], "results", "011-small-qe.txt"}];
records = {"MathematicaVersion=" <> ToString[$Version, InputForm]};

(* r=2: enumerate all five pre-endpoint binary coefficients independently. *)
Clear[b, c]; r = 2; m = 2 r + 2; b[0] = 1;
bits2 = Tuples[{0, 1}, 5];
ans2 = Table[
   Clear[b, c]; b[0] = 1;
   Do[c[k] = bits[[k]], {k, 1, 5}]; c[m] = 1;
   Do[b[k] = Expand[c[k] - s b[k - 1] - If[k >= r, t b[k - r], 0] - If[k >= m, b[k - m], 0]], {k, 1, m}];
   sol = Reduce[0 < s < 1 && 0 < t < 1 && And @@ Table[b[k] >= 0, {k, 1, m}]
       && b[4] == 0 && b[5] == 0 && b[6] == 0, {s, t}, Reals];
   {bits, sol},
   {bits, bits2}];
If[!And @@ (Last[#] === False & /@ ans2), Print["FAILED r=2"]; Exit[2]];
AppendTo[records, "r=2 branches=32 results=" <> ToString[InputForm[ans2]]];

(* r>=3: the endpoint-prefix deductions leave only u=c_(2r). *)
Do[
  Clear[b, c]; m = 2 r + 2; b[0] = 1;
  ans = Table[
    Clear[b, c]; b[0] = 1;
    Do[c[k] = Which[k == 2 r, v, k == m, 1, True, 1], {k, 1, m}];
    Do[b[k] = Expand[c[k] - s b[k - 1] - If[k >= r, t b[k - r], 0] - If[k >= m, b[k - m], 0]], {k, 1, m}];
    sol = Reduce[0 < s < 1 && 0 < t < 1 && And @@ Table[b[k] >= 0, {k, 1, m}]
       && b[r + 2] == 0 && b[2 r + 1] == 0 && b[2 r + 2] == 0,
       {s, t}, Reals];
    {v, sol},
    {v, {0, 1}}];
  If[!And @@ (Last[#] === False & /@ ans), Print["FAILED r=", r]; Exit[2]];
  AppendTo[records, "r=" <> ToString[r] <> " branches=" <> ToString[InputForm[ans]]],
  {r, 3, 8}];
Scan[Print, records];
Export[out, StringRiffle[records, "\n"] <> "\n", "Text"];
Print["RESULT_FILE=", out];
