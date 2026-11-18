from fractions import Fraction as Q

# |z| <= |u|/4 gives Re((u-z)^2)/u^2 >= 1/2.
assert (1-Q(1,4))**2-Q(1,4)**2 == Q(1,2)
# Small disk within the uniform holomorphic disk.
assert Q(1,32)/Q(3,8) == Q(1,12)
print("Both exact disk inequalities verified.")
