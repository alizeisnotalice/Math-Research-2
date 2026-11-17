$AllowInternet=False;
emit[x_]:=(WriteString[$Output,ExportString[x,"RawJSON","Compact"->True],"\n"];Flush[$Output]);
gpu[x_]:=GPUArray[NumericArray[N[x],"Real32"]];
ensure[x_]:=If[!TrueQ[GPUArrayQ[x]],emit[<|"event"->"gpu_fallback_rejected"|>];Exit[3]];
dc=gpu[data["DC"]]; ds=gpu[data["DS"]]; ensure[dc];ensure[ds];
started=AbsoluteTime[];
emit[<|"event"->"started","pid"->$ProcessID,"version"->$Version,"tasks"->Length[data["tasks"]],"batch"->data["batch"]|>];
Do[
 lr=gpu[task["LR"]];li=gpu[task["LI"]];ng=task["groups"]; total=task["masks"]; start=task["resume_start_mask"];
 If[!IntegerQ[start]||start<0||start>=total,emit[<|"event"->"invalid_resume_start","task"->task["id"],"start"->start,"masks"->total|>];Exit[4]];
 emit[<|"event"->"task_resume_v2_started","task"->task["id"],"first"->start,"last_possible"->(total-1)|>];
 Do[
  codes=Range[offset,Min[total-1,offset+data["batch"]-1]];
  mask=Developer`ToPackedArray[N[Prepend[IntegerDigits[#,2,ng-1],0]& /@ codes]];
  mm=gpu[mask];nn=gpu[1.-mask];
  ar=mm.lr;ai=mm.li;br=nn.lr;bi=nn.li;
  Scan[ensure,{ar,ai,br,bi}];
  avR=Exp[ar] Cos[ai];avI=Exp[ar] Sin[ai];bvR=Exp[br] Cos[bi];bvI=Exp[br] Sin[bi];
  ac=avR.dc+avI.ds;bc=bvR.dc+bvI.ds;
  ia=avI.dc-avR.ds;ib=bvI.dc-bvR.ds;Scan[ensure,{ac,bc,ia,ib}];
  aa=Normal[ac];bb=Normal[bc];imA=Normal[ia];imB=Normal[ib];
  scores=MapThread[Max[0.,-Min[#1],-Min[#2],Max[#1]-1.,Max[#2]-1.,Max[Abs[#3]],Max[Abs[#4]]]&,{aa,bb,imA,imB}];
  (* Remove constant/trivial assignments, preserving every other numerical score. *)
  If[offset==0,scores[[1]]=Infinity];
  best=Ordering[scores,Min[4,Length[scores]]];
  emit[<|"event"->"batch","task"->task["id"],"first"->offset,"last"->Last[codes],"count"->Length[codes],"gpu_arrays_verified"->True,
    "cases"->Table[<|"mask"->codes[[i]],"score"->scores[[i]],"A"->aa[[i]],"B"->bb[[i]]|>,{i,best}],"seconds"->AbsoluteTime[]-started|>];
  Clear[mm,nn,ar,ai,br,bi,avR,avI,bvR,bvI,ac,bc,ia,ib,aa,bb,imA,imB,mask];
 ,{offset,start,total-1,data["batch"]}];
 Clear[lr,li];emit[<|"event"->"task_resume_v2_complete","task"->task["id"],"first"->start,"end"->total|>];
,{task,data["tasks"]}];
Clear[dc,ds];emit[<|"event"->"complete","seconds"->AbsoluteTime[]-started,"pid"->$ProcessID|>];Exit[];
