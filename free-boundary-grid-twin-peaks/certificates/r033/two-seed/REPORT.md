# Two-seed completion comparison: exact finite counterexamples

Let B_(p,r) count descending vertex orders with ordered nonadjacent top roots
p,r and with every later vertex adjacent to an earlier vertex.

Full 9! enumeration on the 3-by-3 free grid gives
B_((0,0),(2,1))=2100 at Manhattan distance3, and
B_((0,0),(2,2))=1920 at distance4. The n4 complete subset recurrence
F(V)=1, F(S)=sum_{v outside S, v adjacent S} F(S+v) yields
B_((1,0),(1,3))=2651194000 at distance3 and
B_((1,0),(3,2))=2523561040 at distance4. Both moved roots are in the same D4
edge orbit. The n3 DP sum81984 agrees with the independent permutation sum.

Thus distance-monotonicity fails even with fixed high root and controlled
type of the other root. No averaged transport inequality is refuted here.

For an event E of the roots and B_E>0, the exact change of measure is
P2(E)=(B_E/B) E_(B_E)L/E_B L. A lower bound on B/A2 would be one sufficient
way to use a positive B-law event; it is not indispensable, because a direct
comparison of these conditional lift means can also suffice. No such
asymptotic comparison or fixed positive far-peak fraction is proved here.

The scripts use Python integers, no random inputs. Replay them in an external
directory and compare the printed JSON with the archived results.
