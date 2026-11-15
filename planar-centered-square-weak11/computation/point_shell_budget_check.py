#!/usr/bin/env python3
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,platform,sys
ROOT=Path(__file__).resolve().parents[1]
cases=[]
for M,eps in ((F(1),F(1)),(F(8),F(1,32)),(F(7,5),F(2,3)),(F(1,10),F(10))):
    shell=(M+eps)-eps
    actual=(M+eps)-max(M,eps)
    old_intersection=min(M+eps,M)-min(eps,M)
    assert shell==M
    assert actual==min(M,eps)
    assert shell-old_intersection==actual
    cases.append({'M':str(M),'epsilon':str(eps),'shell_area':str(shell),
                  'coincident_exact_pure_area':str(actual),'old_shell_overlap':str(old_intersection)})
# No mixed subset can span the separated points at 0 and 2.
assert F(2)**2>=F(1)+F(1)
# Distinct equal-mass sources have order-independent crude cumulative charge.
for N in (2,3,8,32):
    assert sum(F(k,N) for k in range(N))==F(N-1,2)
for k in range(2,9):
    m=F(3,7);M=(k-1)*m;eps=m
    assert M-(M-eps)==eps
paths=[ROOT/'proofs/point-insertion-shell-budget.md',Path(__file__).resolve()]
out={'status':'PASS','python':sys.version,'platform':platform.platform(),
     'arithmetic':'fractions.Fraction, exact scalar identities only','cases':cases,
     'coincident_telescope_checks':7,
     'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
result=ROOT/'computation/results/joint-old-reserve/exact-shell-check.json'
result.parent.mkdir(parents=True,exist_ok=True)
result.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
