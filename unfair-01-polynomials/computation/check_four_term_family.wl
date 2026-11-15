Print["Version: ",$Version];
Do[
 beta[0]=1;
 Do[beta[k]=Expand[1-s beta[k-1]-If[k>=r,t beta[k-r],0]],{k,1,r+1}];
 formula=(1-(-s)^(r+2))/(1+s)+(2s-1)t;
 Print["r=",r," identity=",Together[beta[r+1]-formula],
       " nonpositive=",Reduce[0<s<1 && 0<t<1 && beta[r+1]<=0,{s,t},Reals]];
 Clear[beta],
 {r,2,12}];
