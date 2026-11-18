"""Rational verification of the center density witness."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json

x=Q(5,2)
poly=sum((x**j/Q(factorial(j)) for j in range(13)),Q(0))
bound=poly+(x**13/Q(factorial(13)))/(1-x/14)
assert bound<13
assert Q(39,100)**2<Q(7,44)
leakage=Q(1,64)*Q(571,800)*Q(78,25)/13
max_target=Q(250000,249849)
margin=1+leakage-max_target
assert margin>0
data={"exp_5_over_2_upper":str(bound),"leakage_lower":str(leakage),
      "target_upper":str(max_target),"strict_margin":str(margin)}
Path("center-certificate.json").write_text(json.dumps(data,indent=2)+"\n")
