(* Standalone exact symbolic cross-check, not a substitute for halo geometry. *)
ClearAll[b,w,h,L];
t0=AbsoluteTime[];
f[b_]:=Sqrt[1+b]-Sqrt[b]/2;
alpha=Sqrt[3]/2;
checks=<|
  "normalized_minimum_all_nonnegative_b" -> Resolve[ForAll[b,b>=0,2 Sqrt[1+b]>=Sqrt[b]+Sqrt[3]],Reals],
  "unique_minimizer" -> Resolve[ForAll[b,b>=0,Equivalent[f[b]==alpha,b==1/3]],Reals],
  "minimizer_value" -> FullSimplify[f[1/3]==alpha],
  "square_certificate" -> FullSimplify[(1+3b)^2-12b==(1-3b)^2],
  "scaled_minimum" -> FullSimplify[Sqrt[L+L/3]-Sqrt[L/3]/2==Sqrt[3L]/2,Assumptions->L>0],
  "area_difference" -> FullSimplify[(w+Sqrt[L])(h+Sqrt[L])-(Sqrt[3L]-w)(Sqrt[3L]-h)==(1+Sqrt[3])Sqrt[L](w+h)-2L,Assumptions->L>0],
  "fat_domain_central_containment" -> Resolve[ForAll[{w,L},L>0&&3Sqrt[L]/4<w<Sqrt[L],w>=((Sqrt[3]-1)/2)Sqrt[L]],Reals],
  "central_rectangle_positive" -> Resolve[ForAll[{w,h,L},L>0&&0<=w<Sqrt[L]&&0<=h<Sqrt[L],(Sqrt[3L]-w)(Sqrt[3L]-h)>0],Reals],
  "fat_surplus_area_strict_cap" -> Resolve[ForAll[{w,h,L},L>0&&3Sqrt[L]/4<w<Sqrt[L]&&3Sqrt[L]/4<h<Sqrt[L],(1+Sqrt[3])Sqrt[L](w+h)-2L<2Sqrt[3]L],Reals]
|>;
result=<|"schema"->"r068-central-collar-symbolic-v1","status"->If[And@@(TrueQ/@Values[checks]),"PASS","FAIL"],"version"->$Version,"system_id"->$SystemID,"checks"->checks,"exact_minimum"->ToString[alpha,InputForm],"unique_minimizer"->"1/3","algebra_scope"->"Exact scalar minimum, scaling, containment inequality and area identity. No geometric source enumeration, general selection estimate or C2 improvement is asserted.","kernel_wall_seconds"->N[AbsoluteTime[]-t0]|>;
outfile=FileNameJoin[{DirectoryName[$InputFileName],"r068-central-collar-symbolic.json"}];
written=Export[outfile,result,"RawJSON"];
Print[result["status"]];
If[written===$Failed||result["status"]!="PASS",Exit[1],Exit[0]];
