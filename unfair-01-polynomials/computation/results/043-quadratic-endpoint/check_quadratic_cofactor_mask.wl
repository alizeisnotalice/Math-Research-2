(* Exact obstruction to the false inference G>=0 and GP>=0 => P>=0. *)
Clear[x];
alpha=2-Sqrt[2];
P=1-alpha*x+x^2;
G=1+x;
A=Expand[G*P];
Pbar=Expand[P/.Sqrt[2]->-Sqrt[2]];
H=Expand[P*Pbar];
result=<|"P"->P,"G"->G,"A"->A,"Pbar"->Pbar,"normH"->H,
 "P_x1_negative"->TrueQ[Coefficient[P,x,1]<0],
 "A_coefficients_nonnegative"->And@@Thread[CoefficientList[A,x]>=0],
 "A_constant_and_leading_one"->(Coefficient[A,x,0]==1&&Coefficient[A,x,Exponent[A,x]]==1),
 "A_middle_coefficients"->Rest[Most[CoefficientList[A,x]]]|>;
Print[result];
If[!TrueQ[result["P_x1_negative"]&&result["A_coefficients_nonnegative"]&&result["A_constant_and_leading_one"]],Exit[1]];
