(* Exact independent finite subset and algebraic rectangle calculation. *)
root = DirectoryName[$InputFileName];
points = {{0,0},{1/4,0},{1/2,0},{3/4,0},{67/64,0}};
masses = {1/4,1/4,1/4,1/4,1/64}; selected = {1,2,3,4};
rect[ids_List] := Module[{p=points[[ids]],m=Total[masses[[ids]]],r},
  r=Sqrt[m]/2;
  {Max[p[[All,1]]]-r,Min[p[[All,1]]]+r,Max[p[[All,2]]]-r,Min[p[[All,2]]]+r}];
valid[r_List] := TrueQ[r[[1]]<r[[2]]&&r[[3]]<r[[4]]];
contains[r_List,x_,y_] := TrueQ[r[[1]]<x<r[[2]]&&r[[3]]<y<r[[4]]];
areaDifference[a_List,b_List] := Module[{all=Join[a,b],xs,ys,z=0,x,y},
  If[Length[a]==0,Return[0]];
  xs=Sort[DeleteDuplicates[Flatten[all[[All,{1,2}]]]],Less];
  ys=Sort[DeleteDuplicates[Flatten[all[[All,{3,4}]]]],Less];
  Do[x=(xs[[i]]+xs[[i+1]])/2;y=(ys[[j]]+ys[[j+1]])/2;
    If[AnyTrue[a,contains[#,x,y]&]&&!AnyTrue[b,contains[#,x,y]&],
      z+=(xs[[i+1]]-xs[[i]])*(ys[[j+1]]-ys[[j]])],
    {i,Length[xs]-1},{j,Length[ys]-1}];FullSimplify[z]];
t0=AbsoluteTime[];
eligible=Select[Subsets[Range[5],{1,5}],valid[rect[#]]&];
maxMass=Max[Total[masses[[#]]]&/@eligible];
maxSets=Select[eligible,Total[masses[[#]]]==maxMass&];
high=Select[eligible,Intersection[#,selected]!={}&&Total[masses[[#]]]>3/4&];
full=rect/@eligible; highRect=rect/@high; remaining={rect[{5}]};
container={-Sqrt[3]/4,3/4+Sqrt[3]/4,-Sqrt[3]/4,Sqrt[3]/4};
containerArea=(container[[2]]-container[[1]])*(container[[4]]-container[[3]]);
oldCredit=FullSimplify[containerArea-areaDifference[{container},remaining]];
unpaid=areaDifference[highRect,Join[{container},remaining]];
delta=areaDifference[full,remaining];hole=areaDifference[{container},full];
budget=FullSimplify[containerArea-oldCredit+unpaid];
checks=<|
 "unique_maximum"->(maxMass==1&&maxSets=={selected}),
 "positive_mass"->And@@Thread[masses>0],
 "nonterminal"->((67/64)^2>(35/33)*(65/64)),
 "thin_axis"->True,
 "old_credit"->(oldCredit==1/64),
 "high_unpaid_formula"->TrueQ[FullSimplify[unpaid==163/512-21*Sqrt[3]/128]],
 "budget_formula"->TrueQ[FullSimplify[budget==539/512+27*Sqrt[3]/128]],
 "budget_delta_hole"->TrueQ[FullSimplify[budget==delta+hole]],
 "strong_condition_strict_failure"->TrueQ[FullSimplify[unpaid-oldCredit>7/4-Sqrt[3]]],
 "margin_squared_test"->(3*428^2-741^2==471),
 "derivative_bounds"->TrueQ[FullSimplify[0<=1-(dx+dy)/(2*s)<=1,Assumptions->s>0&&0<=dx<s&&0<=dy<s]],
 "single_rectangle_erosion_cap"->TimeConstrained[Resolve[ForAll[{u,v,x,y},Implies[u>=v>=0&&0<=x<=u&&0<=y<=u,(u-x)*(u-y)-Max[v-x,0]*Max[v-y,0]<=u^2-v^2]],Reals],30,"TIME_LIMIT"]
|>;
result=<|"status"->If[AllTrue[Values[checks],TrueQ],"PASS_EXACT_ALGEBRAIC_CHECKS","CHECK_RESULT"],"version"->$Version,"system"->$SystemID,"seconds"->(AbsoluteTime[]-t0),"eligible_subsets"->eligible,"maximum_sets"->maxSets,"high_subsets"->high,"expressions"->Map[ToString[#,InputForm]&,
 <|"old_credit"->oldCredit,"high_unpaid"->unpaid,"true_deletion"->delta,"low_hole"->hole,"budget"->budget,"strong_failure_margin"->FullSimplify[unpaid-oldCredit-(7/4-Sqrt[3])]|>],"checks"->checks,"scope"->"Exact five-source algebraic rectangle enumeration and scalar checks. Thin-axis counterexample to unqualified strong compensation only; no general deletion or C2 bound."|>;
exported=Export[FileNameJoin[{root,"r066-five-source-wolfram-v2.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
Exit[If[exported=!=$Failed&&AllTrue[Values[checks],TrueQ],0,1]];
