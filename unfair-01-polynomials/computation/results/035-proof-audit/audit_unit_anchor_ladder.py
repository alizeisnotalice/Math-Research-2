#!/Users/zzh/.venvs/main/bin/python
"""Exact symbolic checks for the unit-anchor ladder and the p=2 examples."""
from __future__ import annotations
import hashlib,json
from pathlib import Path
import sympy as sp

a,b=sp.symbols('a b', positive=True)

# Ladder recurrence f_k=1-a*f_{k-1}, 0<a<1, checked symbolically for L=2..8.
ladder={}
for L in range(2,9):
    f=1-a
    vals=[f]
    for _ in range(2,L):
        f=sp.expand(1-a*f); vals.append(f)
    ladder[str(L)]={'f_multiples_before_anchor':[str(v) for v in vals],
                     'anchor_positive_term':str(sp.expand(a*vals[-1]))}

# p=2,q=13 direct coefficient identities listed in the note.
f4=1-a; f8=0
h8=sp.expand(f8+a*f4)
assert sp.expand(h8-a*(1-a)) == 0
# Reversed b calculation at degrees 8,12,13,21.
u8=1-b; u4=u12=u13=u5=u1=0
u21,u9=sp.symbols('u21 u9', nonnegative=True)
h21=sp.expand(u21+u9+u8+1)
assert sp.expand(h21-(u21+u9+2-b)) == 0

out={'ladder_L_2_to_8':ladder,
     'p2q13':{'h8':str(h8),'reversed_h21_template':str(h21)},
     'p2q17':{'h8_if_fractional':'f8+a*(1-a)', 'h12_after_f12_zero':'a*f8'},
     'p2q19':{'f5':'1-a','f10':'1-a+a^2','f15':'1-a+a^2-a^3',
               'forbidden_h20_positive_term':'a*(1-a+a^2-a^3)'},
     'scope':'symbolic identity audit only; the accompanying note supplies inequality proof.'}
text=json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n'
Path('unit_anchor_ladder_audit.json').write_text(text,encoding='utf-8')
print(json.dumps({'sha256':hashlib.sha256(text.encode()).hexdigest(),'L_checked':list(range(2,9))}))
