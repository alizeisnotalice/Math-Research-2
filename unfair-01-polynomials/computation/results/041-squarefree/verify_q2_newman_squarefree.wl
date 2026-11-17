(* Independent exact verifier for the theorem:
   no K in {0,1}[x] with K(0)=1 is divisible by (1+x^2+x^3)^2.
   It imports neither the numerical pilot nor its checkpoint. *)
Clear[x,y,step,poly,key];
q=1+x^2+x^3;
states={{1,0,0,0,0,0},{0,1,0,0,0,0},{1,1,0,0,0,0},
 {0,0,1,0,0,0},{1,0,1,0,0,0},{0,1,1,0,0,0},{1,1,1,0,0,0},
 {0,0,0,1,0,0},{1,0,0,1,0,0},{0,0,1,1,0,0},{1,0,1,1,0,0},
 {0,1,1,1,0,0},{1,1,1,1,0,0},{0,0,0,1,1,0},{1,0,0,1,1,0},
 {0,1,0,1,1,0},{1,1,0,1,1,0}};
step[a_List,bit_Integer]:=With[{t=a[[6]]},{bit-t,a[[1]],a[[2]]-2t,a[[3]]-2t,a[[4]]-t,a[[5]]-2t}];
poly[a_List]:=Sum[a[[i+1]]x^i,{i,0,5}];
key[a_List]:=StringRiffle[ToString/@a,","];
seen=AssociationThread[key/@states->Range[Length[states]]];
region=(q/.x->y)==0 && -22/15<y<-10/7;
rootAudit=Resolve[Exists[y,region],Reals] &&
  Resolve[ForAll[y,region\[Implies] D[q/.x->y,y]>0],Reals];
edges=Flatten[Table[{a,b,step[a,b]},{a,states},{b,0,1}],1];
transitionAudit=And@@Flatten[Table[
 PolynomialRemainder[x poly[a]+b,q^2,x]===poly[step[a,b]],{a,states},{b,0,1}]];
outside=Select[edges,!KeyExistsQ[seen,key[#[[3]]]]&];
cert[edge_]:=Module[{z=edge[[3]],v,dv},
 v=PolynomialRemainder[poly[z],q,x]/.x->y;
 dv=PolynomialRemainder[D[poly[z],x],q,x]/.x->y;
 Resolve[ForAll[y,region\[Implies](v^2>9 || dv^2>36)],Reals]
];
edgeAudit=cert/@outside;
result=<|"rootIntervalAndMonotonicity"->rootAudit,"stateCount"->Length[states],
 "zeroNotInStateSet"->!KeyExistsQ[seen,key[{0,0,0,0,0,0}]],
 "transitionFormulaModuloQSquared"->transitionAudit,"edgeCount"->Length[edges],"outsideEdgeCount"->Length[outside],
 "allOutsideEdgesViolateValueOrDerivativeBound"->And@@edgeAudit,
 "outsideEdges"->outside,"edgeAudit"->edgeAudit|>;
Print[result];
If[!TrueQ[rootAudit && transitionAudit && result["zeroNotInStateSet"] && result["allOutsideEdgesViolateValueOrDerivativeBound"]],Exit[1]];
