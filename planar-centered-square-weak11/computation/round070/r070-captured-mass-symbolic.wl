(* Exact scalar checks only; the source geometry is proved separately. *)
ClearAll[a,m,u,v,r,c,d]; t0=AbsoluteTime[];
f[a_,m_]:=Sqrt[1+m-a]-Sqrt[m]/2;
expr=D[f[a,m],m];
checks=<|
 "derivative_rationalized"->FullSimplify[expr==(3m+a-1)/(4 Sqrt[m] Sqrt[1+m-a](2 Sqrt[m]+Sqrt[1+m-a])),Assumptions->0<a<1&&a<=m<=1],
 "square_lower_bound_identity"->Expand[(v-u/2)^2-3(v^2-u^2)/4-(v/2-u)^2]===0,
 "small_a_minimizer_feasible"->Resolve[ForAll[a,0<a<=1/4,a<=(1-a)/3<=1],Reals],
 "small_a_minimum_value"->FullSimplify[f[a,(1-a)/3]==Sqrt[3(1-a)]/2,Assumptions->0<a<=1/4],
 "large_a_derivative_nonnegative"->Resolve[ForAll[{a,m},1/4<=a<1&&a<=m<=1,3m+a-1>=0],Reals],
 "large_a_minimum_value"->FullSimplify[f[a,a]==1-Sqrt[a]/2,Assumptions->1/4<=a<1],
 "branch_values_match"->FullSimplify[Sqrt[3(1-1/4)]/2==1-Sqrt[1/4]/2],
 "mass_radius_equivalence"->Resolve[ForAll[{a,r,c},0<a<1&&r>0&&c>=0,Equivalent[c>Sqrt[1+4r^2-a]-r,a>1+3r^2-2c r-c^2]],Reals],
 "lower_c_inverse"->Resolve[ForAll[{a,c},1/4<=a<1&&1/2<c<=3/4,Equivalent[c>1-Sqrt[a]/2,a>4(1-c)^2]],Reals],
 "upper_c_inverse"->Resolve[ForAll[{a,c},0<a<=1/4&&3/4<=c<=Sqrt[3]/2,Equivalent[c>Sqrt[3(1-a)]/2,a>1-4c^2/3]],Reals],
 "small_a_cannot_lower_c"->Resolve[ForAll[{a,c},0<a<=1/4&&1/2<c<=3/4,c<=Sqrt[3(1-a)]/2],Reals],
 "large_a_always_upper_c"->Resolve[ForAll[{a,c},1/4<a<1&&3/4<=c<=Sqrt[3]/2,c>1-Sqrt[a]/2],Reals],
 "clamp_endpoint_with_distance"->Resolve[ForAll[{a,d},0<a<1&&0<=d<1/2,Max[a,4d^2,(1-a)/3]<=1],Reals],
 "no_capture_at_half_core"->Resolve[ForAll[a,0<a<1,If[a<=1/4,Sqrt[3(1-a)]/2,1-Sqrt[a]/2]>1/2],Reals]
|>;
res=<|"schema"->"r070-captured-mass-symbolic-v1","version"->$Version,"system_id"->$SystemID,"checks"->checks,"status"->If[And@@(TrueQ/@Values[checks]),"PASS","FAIL"],"kernel_wall_seconds"->N[AbsoluteTime[]-t0],"scope"->"Normalized L=1 scalar minimum, derivative, radical-to-polynomial equivalence and inverse branches. Analytic source geometry and area bounds require separate proof. No general C2 bound asserted."|>;
Export[FileNameJoin[{DirectoryName[$InputFileName],"r070-captured-mass-symbolic.json"}],res,"RawJSON"];
Print[res["status"]]; If[res["status"]==="PASS",Exit[0],Exit[1]];
