(* GPU-resident dense coefficient optimization.  Default is a small smoke;
   set runMode="production" only after its GPUArrayQ receipts are inspected. *)
$AllowInternet=False;
emit[x_]:=(WriteString[$Output,ExportString[x,"RawJSON","Compact"->True],"\n"];Flush[$Output]);
gpu[x_]:=GPUArray[NumericArray[N[x],"Real32"]];
ensure[x_]:=If[!TrueQ[GPUArrayQ[x]],emit[<|"event"->"gpu_fallback_rejected"|>];Exit[3]];
runMode="production";m=60;n=100;l=256;seed=2026111617;reportEvery=100;maxSeconds=4*60*60;
{starts,steps}=If[runMode=="smoke",{128,4},{8192,12000}];
lr=.003;b1=.9;b2=.999;eps=10.^-7;clen=m+n+1;
If[clen>l,emit[<|"event"->"bad_dft_length"|>];Exit[2]];
fr[q_]:=N[Table[Cos[2. Pi r s/l],{r,0,q-1},{s,0,l-1}]];
fi[q_]:=N[Table[-Sin[2. Pi r s/l],{r,0,q-1},{s,0,l-1}]];
ir[q_]:=N[Table[Cos[2. Pi r s/l]/l,{r,0,l-1},{s,0,q-1}]];
ii[q_]:=N[Table[-Sin[2. Pi r s/l]/l,{r,0,l-1},{s,0,q-1}]];
{faR,faI,fbR,fbI,fcR,fcI,iaR,iaI,ibR,ibI,icR,icI}=gpu /@ {fr[m+1],fi[m+1],fr[n+1],fi[n+1],fr[clen],fi[clen],ir[m+1],ii[m+1],ir[n+1],ii[n+1],ir[clen],ii[clen]};
ones=gpu[ConstantArray[1.,clen]];Scan[ensure,{faR,faI,fbR,fbI,fcR,fcI,iaR,iaI,ibR,ibI,icR,icI,ones}];
fixed[len_]:=ReplacePart[ConstantArray[0.,len],{1->1.,2->.5,len->1.}];
free[len_]:=ReplacePart[ConstantArray[1.,len],{1->0.,2->0.,len->0.}];
fixA=gpu[ConstantArray[fixed[m+1],starts]];fixB=gpu[ConstantArray[fixed[n+1],starts]];
freeA=gpu[ConstantArray[free[m+1],starts]];freeB=gpu[ConstantArray[free[n+1],starts]];
SeedRandom[seed,Method->"MersenneTwister"];za=gpu[RandomReal[{-2.,2.},{starts,m+1}]];zb=gpu[RandomReal[{-2.,2.},{starts,n+1}]];
ma=0. za;mb=0. zb;va=0. za;vb=0. zb;Scan[ensure,{fixA,fixB,freeA,freeB,za,zb,ma,mb,va,vb}];
sig[z_]:=1./(1.+Exp[-z]); started=AbsoluteTime[];timeout=False;completed=0;
emit[<|"event"->"started","pid"->$ProcessID,"mode"->runMode,"m"->m,"n"->n,"dft_length"->l,"starts"->starts,"steps_requested"->steps,"seed"->seed,"scope"->"finite nonconvex optimization, not exhaustive search/proof/certificate"|>];
Do[
 If[AbsoluteTime[]-started>=maxSeconds,timeout=True;Break[]];
 a=fixA+freeA sig[za];b=fixB+freeB sig[zb];
 ar=Dot[a,faR];ai=Dot[a,faI];br=Dot[b,fbR];bi=Dot[b,fbI];cr=ar br-ai bi;ci=ar bi+ai br;c=Dot[cr,icR]+Dot[ci,icI];
 loss=Dot[c^2 (c-1.)^2,ones];e=2. c (c-1.) (2. c-1.);er=Dot[e,fcR];ei=Dot[e,fcI];
 (* IFFT(E Conjugate[B]) gives Sum_j E_(r+j) B_j.  For r<=m,
    r+j<=m+n<l, so the selected correlation coordinates do not wrap. *)
 grA=er br+ei bi;giA=ei br-er bi;gA=Dot[grA,iaR]+Dot[giA,iaI];
 grB=er ar+ei ai;giB=ei ar-er ai;gB=Dot[grB,ibR]+Dot[giB,ibI];
 sa=sig[za];sb=sig[zb];gzA=gA freeA sa (1.-sa);gzB=gB freeB sb (1.-sb);
 ma=b1 ma+(1.-b1) gzA;mb=b1 mb+(1.-b1) gzB;va=b2 va+(1.-b2) gzA^2;vb=b2 vb+(1.-b2) gzB^2;
 mhA=ma/(1.-b1^step);mhB=mb/(1.-b1^step);vhA=va/(1.-b2^step);vhB=vb/(1.-b2^step);
 denA=Sqrt[vhA]+eps;denB=Sqrt[vhB]+eps;za=za-lr mhA/denA;zb=zb-lr mhB/denB;
 Scan[ensure,{a,b,ar,ai,br,bi,cr,ci,c,loss,e,er,ei,grA,giA,gA,grB,giB,gB,sa,sb,gzA,gzB,ma,mb,va,vb,mhA,mhB,vhA,vhB,denA,denB,za,zb}];
 If[Mod[step,reportEvery]==0,scores=Normal[loss];best=Ordering[scores,Min[4,Length[scores]]];emit[<|"event"->"progress","step"->step,"score_count"->Length[scores],"minimum_loss"->Min[scores],"best_starts"->(best-1),"seconds"->AbsoluteTime[]-started|>];Clear[scores,best]];
 completed=step;
 Clear[a,b,ar,ai,br,bi,cr,ci,c,loss,e,er,ei,grA,giA,gA,grB,giB,gB,sa,sb,gzA,gzB,mhA,mhB,vhA,vhB,denA,denB];
,{step,1,steps}];
(* Terminal evaluation; score vector plus four selected coefficient rows are the only CPU transfers. *)
a=fixA+freeA sig[za];b=fixB+freeB sig[zb];ar=Dot[a,faR];ai=Dot[a,faI];br=Dot[b,fbR];bi=Dot[b,fbI];cr=ar br-ai bi;ci=ar bi+ai br;c=Dot[cr,icR]+Dot[ci,icI];loss=Dot[c^2 (c-1.)^2,ones];Scan[ensure,{a,b,ar,ai,br,bi,cr,ci,c,loss}];
scores=Normal[loss];best=Ordering[scores,Min[4,Length[scores]]];sel=SparseArray[(Rule[#,1.]& /@ Transpose[{Range[Length[best]],best}]),{Length[best],starts}];selGPU=gpu[Normal[sel]];pa=selGPU.a;pb=selGPU.b;pc=selGPU.c;Scan[ensure,{selGPU,pa,pb,pc}];outA=Normal[pa];outB=Normal[pb];outC=Normal[pc];
emit[<|"event"->"final","steps_completed"->completed,"timed_out"->timeout,"score_count"->Length[scores],"cases"->Table[<|"start"->best[[i]]-1,"loss"->scores[[best[[i]]]],"A"->outA[[i]],"B"->outB[[i]],"C"->outC[[i]]|>,{i,Length[best]}],"seconds"->AbsoluteTime[]-started|>];
Clear[faR,faI,fbR,fbI,fcR,fcI,iaR,iaI,ibR,ibI,icR,icI,ones,fixA,fixB,freeA,freeB,za,zb,ma,mb,va,vb,a,b,ar,ai,br,bi,cr,ci,c,loss,scores,best,sel,selGPU,pa,pb,pc,outA,outB,outC];emit[<|"event"->"complete","timed_out"->timeout,"seconds"->AbsoluteTime[]-started,"pid"->$ProcessID|>];Exit[];
