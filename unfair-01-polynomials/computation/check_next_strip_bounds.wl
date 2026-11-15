(* Independent exact check: all-one endpoint residual and tail identity. *)
Clear[s,y,r];
d=1-2s+3s^2;
t=(1+s^3 y)/((1+s)d);
u2=(1-s^2 y^2)/(1+s);
w=(1-((r+3)s^2+(r+2)s^3)y)/(1+s)^2;
f=Together[u2-t w+t^2(1-3s)];
p=9s^6+6s^5-5s^4+4s^3-3s^2-2s+1;
Print["version=",$Version];
Print["limiting_residual_check=",FullSimplify[(f/.y->0)==p/((1+s)^3 d^2)]];
del=Factor[f-(f/.y->0)];
Print["tail=",InputForm[del]];
aa0=3s^3+s^2-s+1;
cc0=3-3s+3s^2-3s^3;
aa1=3s^6+s^5-s^4+s^3;
cc1=-6s^6-3s^5+2s^4-s^3-3s^2+2s-1;
Print["full_tail_identity_check=",FullSimplify[
 del==s^2*y*(r*aa0+cc0+y*(r*aa1+cc1))/((1+s)^3*d^2)
 ]];
lhs=(1-(-s)^(r+5))/(1+s)-((1-(-s)^(r+3))/((1+s)*d))*Sum[(j+1)(-s)^j,{j,0,4}];
rhs=(-s^3 (5s-4)+(-s)^(r+3)(2s^4-2s^3+2s^2-2s+1))/((1+s)*d);
Print["b_r4_identity=",FullSimplify[lhs==rhs]];
Do[
 Clear[b,c]; b[0]=1;
 Do[c[k]=1,{k,1,2 rr+1}];
 trr=t /. y->(-s)^rr;
 Do[b[k]=Together[c[k]-s b[k-1]-If[k>=rr,trr*b[k-rr],0]],{k,1,2 rr+1}];
 check=FullSimplify[b[2 rr+1]==(u2-t w+t^2(1-3s))/.{r->rr,y->(-s)^rr}];
 Print["recurrence_check_r=",rr," : ",check],
 {rr,3,8}];
