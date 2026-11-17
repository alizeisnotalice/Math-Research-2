(*
  Self-contained exact verifier for the degree-three common-factor result.
  No imports, checkpoints, absolute paths, floating roots, or output files.
  Run: wolframscript -file verify_cubic_degree3_portable.wl
*)
Clear[x,y,require,qOf,candidateQ,runBFS];
require[cond_,label_]:=If[!TrueQ[cond],Print[<|"status"->"FAIL","check"->label|>];Exit[1]];
qOf[{a_,b_}]:=1+a*x+b*x^2+x^3;

(* Necessary envelope for a monic irreducible cubic Q|K, K binary, K(0)=1.
   CountRoots is exact.  With one real root in (-2,-1/2), the complex pair
   has squared modulus -1/r; with three real roots the interval test is direct. *)
candidateQ[{a_,b_}]:=Module[{q=qOf[{a,b}],irr,pos,real,inband},
 irr=Factor[q]===q;
 pos=CountRoots[q,{x,0,Infinity}];
 real=CountRoots[q,{x,-Infinity,Infinity}];
 inband=If[real==1,CountRoots[q,{x,-2,-1/2}]==1,
          If[real==3,CountRoots[q,{x,-2,-1/2}]==3,False]];
 irr && pos==0 && inband
];
allPairs=Flatten[Table[{a,b},{a,-5,5},{b,-5,5}],1];
found=Select[allPairs,candidateQ];
expected=Sort[{{-1,0},{-1,1},{0,-1},{0,1},{1,-1},{1,0},{1,2},{2,1}}];
require[Sort[found]===expected,"eight exact cubic candidates"];
require[And@@Table[CountRoots[qOf[p],{x,-Infinity,Infinity}]==1,{p,found}],"each candidate has one real root"];
require[And@@Table[MemberQ[found,Reverse[p]],{p,found}],"reciprocal pairing"];
require[And@@Table[Expand[x^3(qOf[p]/.x->1/x)]===qOf[Reverse[p]],{p,found}],"reciprocal polynomial identity"];

(* Exact Horner BFS in Z[x]/Q^2.  Any accepted word has a negative real root
   rho in the displayed interval.  Backward recursion gives
   |r(rho)| <= 1/(|rho|-1), |r'(rho)| <= 1/(|rho|-1)^2.
   tailAudit below proves the supplied rational squared bounds dominate these.
*)
runBFS[a_,b_,left_,right_,vb2_,db2_,tag_]:=Module[
 {q=qOf[{a,b}],coef,region,cap=5000,states,seen,n=1,head=1,step,poly,key,keep,cache=<||>,z,kz,
  rootAudit,tailAudit,transitionAudit,result},
 coef=CoefficientList[Expand[q^2],x];
 region=(q/.x->y)==0 && left<y<right;
 rootAudit=TrueQ[Resolve[Exists[y,region],Reals]] &&
   TrueQ[Resolve[ForAll[y,region\[Implies]D[(q/.x->y),y]>0],Reals]];
 tailAudit=TrueQ[Resolve[ForAll[y,region\[Implies](-y>1 && 1/(-y-1)^2<vb2 && 1/(-y-1)^4<db2)],Reals]];
 require[rootAudit,"unique isolated negative root: "<>tag];
 require[tailAudit,"tail bounds dominate value/derivative bounds: "<>tag];
 step[s_,bit_]:=With[{t=s[[6]]},{bit-t*coef[[1]],s[[1]]-t*coef[[2]],s[[2]]-t*coef[[3]],
   s[[3]]-t*coef[[4]],s[[4]]-t*coef[[5]],s[[5]]-t*coef[[6]]}];
 poly[s_]:=Sum[s[[i+1]]x^i,{i,0,5}];
 key[s_]:=StringRiffle[ToString/@s,","];
 keep[s_]:=Module[{kk=key[s],v,dv,ans,res},If[KeyExistsQ[cache,kk],Return[cache[kk]]];
   v=PolynomialRemainder[poly[s],q,x]/.x->y;
   dv=PolynomialRemainder[D[poly[s],x],q,x]/.x->y;
   res=Resolve[Exists[y,region && v^2<=vb2 && dv^2<=db2],Reals];
   require[res===True || res===False,"Resolve returned a non-Boolean state test: "<>tag];
   ans=res;
   cache[kk]=ans;ans
 ];
 states={{1,0,0,0,0,0}};seen=<|key[First[states]]->1|>;
 While[head<=n && n<cap,
   Do[z=step[states[[head]],bit];If[keep[z],kz=key[z];If[!KeyExistsQ[seen,kz],AppendTo[states,z];n++;seen[kz]=n]],{bit,0,1}];
   head++
 ];
 transitionAudit=And@@Flatten[Table[PolynomialRemainder[x poly[s]+bit,q^2,x]===poly[step[s,bit]],{s,states},{bit,0,1}]];
 result=<|"tag"->tag,"rootAudit"->rootAudit,"tailAudit"->tailAudit,"transitionAudit"->transitionAudit,
   "visited"->n,"queueExhausted"->(head>n),"zeroReachable"->KeyExistsQ[seen,key[{0,0,0,0,0,0}]]|>;
 require[transitionAudit,"integer Q^2 transition: "<>tag];
 require[head>n && n<cap,"finite exact BFS exhausted: "<>tag];
 require[!result["zeroReachable"],"zero state excluded: "<>tag];
 result
];

(* One representative from each reciprocal pair. *)
representatives={{-1,0,-4/3,-29/22,256/25,100,"1-y+y3"},
 {-1,1,-2,-11/6,25/16,9/4,"1-y+y2+y3"},
 {0,1,-3/2,-10/7,9,36,"1+y2+y3"},
 {1,2,-2,-7/4,16/9,4,"1+y+2y2+y3"}};
require[And@@Table[MemberQ[found,Take[r,2]],{r,representatives}],"BFS representatives in candidate list"];
reports=runBFS@@@representatives;
Print[<|"status"->"PASS","candidateCount"->Length[found],"candidates"->found,"reciprocalPairsVerified"->True,"representativeReports"->reports|>];
