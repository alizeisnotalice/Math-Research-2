"""Final exact arithmetic for the mixed maximum first-derivative bound."""
from fractions import Fraction

h = Fraction(1, 2**30)
right_positive_lower = Fraction(1, 3**13)
left_negative_upper = 17 * h
margin = right_positive_lower - left_negative_upper
assert margin > 0
print(f"h={h}, positive={right_positive_lower}, negative={left_negative_upper}, margin={margin}")
