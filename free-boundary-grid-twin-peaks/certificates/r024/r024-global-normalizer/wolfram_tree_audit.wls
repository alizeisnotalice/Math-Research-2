(* Independent exact rational audit of the (height 2,height 4) binary bridge. *)
ClearAll[recPeak, recNoPeak, trunc, z, x];
trunc[e_] := Normal[Series[Expand[e], {z, 0, 2}]];
binaryEdges[offset_, h_] := Module[{m = 2^(h + 1) - 1},
  Flatten[Table[With[{a = offset + i},
    ({a, offset + #} & /@ Select[{2 i + 1, 2 i + 2}, # < m &])],
    {i, 0, m - 1}], 1]];
edges = Join[binaryEdges[0, 2], binaryEdges[7, 4], {{0, 7}}];
n = 38;
adj = Table[{}, {n}];
Scan[(AppendTo[adj[[#[[1]] + 1]], #[[2]]]; AppendTo[adj[[#[[2]] + 1]], #[[1]]]) &, edges];
recPeak[v_, parent_, xx_] := Module[{aa = 1, bb = 1, ss, child, rw, aw, bw},
  Do[If[w != parent,
    ss = Unique["s"]; rw = recPeak[w, v, ss]; aw = rw[[1]]; bw = rw[[2]];
    child = Integrate[aw, {ss, 0, 1}] + (z - 1) Integrate[bw, {ss, xx, 1}];
    aa = trunc[aa child]; bb = trunc[bb Integrate[aw, {ss, 0, xx}]];
  ], {w, adj[[v + 1]]}]; {aa, bb}
];
rootPair = recPeak[0, -1, x];
peakPoly = trunc[rootPair[[1]] + (z - 1) rootPair[[2]]];
a2 = Factorial[n] Coefficient[Integrate[peakPoly, {x, 0, 1}], z, 2];
componentNoPeak[vertices_, marked_] := Module[{root = First[marked], walk},
  walk[v_, parent_, xx_] := Module[{aa = 1, bb = 1, ss, rw, aw, bw, child},
    Do[If[w != parent && MemberQ[vertices, w],
      ss = Unique["s"]; rw = walk[w, v, ss]; aw = rw[[1]]; bw = rw[[2]];
      child = If[MemberQ[marked, w], Integrate[aw, {ss, 0, 1}],
        Integrate[aw, {ss, 0, 1}] - Integrate[bw, {ss, xx, 1}]];
      aa = Expand[aa child]; bb = Expand[bb Integrate[aw, {ss, 0, xx}]];
    ], {w, adj[[v + 1]]}]; {aa, bb}
  ]; Integrate[walk[root, -1, x][[1]], {x, 0, 1}]
];
pairQ[u_, v_] := Module[{rem, comps = {}, seed, stack, c, w, q = 1, marked},
  rem = Complement[Range[0, n - 1], {u, v}];
  While[Length[rem] > 0,
    seed = First[rem]; stack = {seed}; c = {seed}; rem = DeleteCases[rem, seed];
    While[Length[stack] > 0,
      w = Last[stack]; stack = Most[stack];
      Do[If[MemberQ[rem, qv], rem = DeleteCases[rem, qv]; AppendTo[c, qv]; AppendTo[stack, qv]],
        {qv, adj[[w + 1]]}]
    ]; AppendTo[comps, c]
  ];
  Do[marked = Select[c, MemberQ[adj[[# + 1]], u] || MemberQ[adj[[# + 1]], v] &];
    q *= componentNoPeak[c, marked], {c, comps}]; q
];
bp = Sum[If[u != v && !MemberQ[adj[[u + 1]], v], pairQ[u, v], 0], {u, 0, n - 1}, {v, 0, n - 1}]/(n (n - 1));
b2 = Factorial[n] bp;
Print[<|"vertices" -> n, "A2" -> a2, "Btop2" -> b2, "ratio" -> Together[a2/b2], "ratio_float" -> N[a2/b2, 20]|>];
