(* Scalar identities only. Geometric capture and all-square interpretation
   are analytic obligations in the R73 proof, not supplied by this file. *)
ClearAll[n,k,l,d,r,z];
t=n/(n+1);c=t/(1+2d);
ak=(k(n+1)-n)/n^2;al=(l(n+1)-n)/n^2;
g=(n-k)(n-l)/n^3;
excess=k l/n^2+c((r+2d)^2-(ak+2d)(al+2d))-r^2;
v=(1+1/n^2)^2+t(3-6/n^2);
domain=n>=2&&1<=k<=n&&1<=l<=n&&0<=d<=1/(2n)&&r>=1;
checks=<|
 "unit_excess_factor"->(Together[(excess/.r->1)+(1-2n d)g/(1+2d)]===0),
 "general_radius_factor"->(Together[excess-(excess/.r->1)+(r-1)((1-c)(r+1)-4c d)]===0),
 "radius_bracket_factor"->(Together[(1-c)(r+1)-4c d-((1-c)(r-1)+2(1-t))]===0),
 "zero_capture_cap"->FullSimplify[c(1+2d)^2<=1,Assumptions->n>=2&&0<=d<=1/(2n)],
 "unit_excess_nonpositive"->FullSimplify[-(1-2n d)g/(1+2d)<=0,Assumptions->domain],
 "radius_bracket_positive"->FullSimplify[(1-c)(r-1)+2(1-t)>0,Assumptions->domain],
 "density_between_zero_and_one"->FullSimplify[0<c<1,Assumptions->domain],
 "corner_threshold_identity"->(Together[4(1/2-1/n^2)^2-((2-1/n^2)^2-7/2)-(1/2+3/n^4)]===0),
 "local_full_halo_limit_above_7half"->FullSimplify[v>7/2,Assumptions->n>=7]
|>;
poly=Expand[2n^4(n+1)Together[v-7/2]]//Together//Expand;
shift=CoefficientList[Expand[poly/.n->z+7],z];
out=<|"scope"->"Exact scalar identities and real-domain inequalities only; not a geometric cap proof or full-halo computation.",
 "version"->$Version,"kernel_pid"->$ProcessID,"checks"->checks,
 "all_true"->And@@Values[checks],"local_full_gap_n7"->ToString[Together[(v-7/2)/.n->7],InputForm],
 "local_full_gap_positive_shift_coefficients"->shift|>;
Export[$ScriptCommandLine[[2]],out,"RawJSON"];
If[TrueQ[out["all_true"]]&&And@@Thread[shift>0],Exit[0],Exit[2]];
