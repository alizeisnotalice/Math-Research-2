ClearAll[x];
A = 1 + 9 x/10 + 181 x^3/1000 + 2421 x^4/40000;
B = 1 + x/10 + 91 x^2/100 + 7371 x^4/8000;
prod = Expand[A B];
checks = <|"coefficients_0_to_5" -> (Coefficient[prod, x, #] & /@ Range[0, 5]),
  "coefficient_6" -> Coefficient[prod, x, 6],
  "endpoint_a5_b5" -> {Coefficient[A, x, 5], Coefficient[B, x, 5]}|>;
If[checks["coefficients_0_to_5"] =!= ConstantArray[1, 6] ||
   !(0 < checks["coefficient_6"] < 1), Exit[1]];
Print[checks];
