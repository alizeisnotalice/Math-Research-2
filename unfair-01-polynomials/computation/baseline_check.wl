(* Exact tool smoke and foundational finite sanity checks; not a general proof. *)
Print["Version: ", $Version];
Print["ComplementAndZero: ", Reduce[a+b==1 && a b==0 && 0<a<1 && 0<b<1,{a,b},Reals]];
Print["SymmetricSupportObstruction: ", Reduce[0<a<1 && 0<b<=1 && 1+a b<=1,{a,b},Reals]];
Print["QuadraticFactors: ", Reduce[0<=a<=1 && 0<=b<=1 && (0<a<1 || 0<b<1) && (a+b==0 || a+b==1) && (2+a b==0 || 2+a b==1),{a,b},Reals]];
Print["HarePrintedTotalDifference: ", 7130374+7438+2900+975-7141686];
