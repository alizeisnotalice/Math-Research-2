(* Exact scalar identities only; geometric and measure arguments are separate. *)
Clear[n,k,l,z];
t=n/(n+1); h=(n+1)/n^2; eps=1/n^2;
a[k_]=(k-1)h+eps;
f[k_,l_]=k l/n^2-t a[k]a[l];
cArea=(2-eps)^2; jArea=(1+eps)^2;
bound=jArea+(1-t)(cArea-jArea);
checks=<|
 "cap_gain_increment"->TrueQ[FullSimplify[f[k+1,l]-f[k,l]-(n-l)/n^3]==0],
 "cap_gain_full_grid"->TrueQ[FullSimplify[f[n,n]-(1-t)]==0],
 "sharp_constant_density_unit_window"->TrueQ[FullSimplify[(n-1)/n+t h-1]==0],
 "outside_collar_area"->TrueQ[FullSimplify[cArea-jArea-(3-6/n^2)]==0],
 "atomic_mass_error_coefficient"->TrueQ[FullSimplify[t(12z+4z^2)/(1-t)-n(12z+4z^2)]==0],
 "small_side_mass_bound"->TrueQ[FullSimplify[(1+2z)^2<=1+n(12z+4z^2),Assumptions->n>=2&&z>=0]],
 "atomic_deletion_limit"->TrueQ[Limit[bound,n->Infinity]==1],
 "n5_atomic_upper_below_8over5"->TrueQ[(bound/.n->5)<8/5],
 "n17_atomic_upper_below_6over5"->TrueQ[(bound/.n->17)<6/5]
|>;
out=<|"schema"->"r072-diffuse-algebra-v1","version"->$Version,
 "process_id"->$ProcessID,"checks"->checks,
 "n5_atomic_deletion_limsup_bound"->ToString[bound/.n->5,InputForm],
 "n17_atomic_deletion_limsup_bound"->ToString[bound/.n->17,InputForm],
 "scope"->"Exact algebra and one elementary sign inequality; not a proof of the square cap, halo area, or atomic-limit geometry."|>;
Print[ExportString[out,"RawJSON"]];
If[And@@Values[checks],Exit[0],Exit[1]];
