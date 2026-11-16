ClearAll[a,x];
cases={{8, {0,8,16,20,21,24,28,29,32,36,37,40,41}},{12, {0,12,20,21,24,32,33,36,41}}};
results=Table[
 s=cas[[1]];h=Total[x^cas[[2]]];g=(1+x^21)(1+a*x^s+x^20);
 f=Normal[Series[h/g,{x,0,41}]];
 aa=Expand[(f/.x->x^2)(g/.x->x^21)];bb=Expand[(f/.x->x^21)(g/.x->x^2)];
 fv=CoefficientList[f,x];av=Table[Coefficient[aa,x,k],{k,0,83}];bv=Table[Coefficient[bb,x,k],{k,0,881}];
 constraints=And@@Join[Thread[fv>=0],Thread[av>=0],Thread[av<=1],Thread[bv>=0],Thread[bv<=1]];
 timed=AbsoluteTiming[TimeConstrained[Reduce[0<a<1 && constraints,a,Reals],90,$Aborted]];
 <|"s"->s,"HSupport"->cas[[2]],"prefixIdentity"->(Expand[Normal[Series[f*g-h,{x,0,41}]]]===0),"parameterRegion"->ToString[timed[[2]],InputForm],"seconds"->timed[[1]],"FNonzeroCoefficients"->Table[If[Coefficient[f,x,k]=!=0,{k,ToString[Coefficient[f,x,k],InputForm]},Nothing],{k,0,41}]|>,{cas,cases}];
out=<|"version"->$Version,"pid"->$ProcessID,"scope"->"Only prefixes F,H through degree41, A below84 and B below882; not a complete factorization", "results"->results|>;
Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram-prefix-audit.json"}],out,"RawJSON"];Print[InputForm[out]];Exit[];
