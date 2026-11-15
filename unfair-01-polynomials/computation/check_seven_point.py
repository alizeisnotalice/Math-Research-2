#!/usr/bin/env python3
import json
import sympy as S

x=S.Symbol('x')
cases=[([0,1,2,4,7],[0,1,2,4,8,10,13],[4,5]),
       ([0,1,2,3,5,11],[0,1,3,4,5,7,13],[15,18,9]),
       ([0,1,2,4,6,11],[0,1,3,4,6,8,13],[1,19,9])]
outputs=[]
for index,(I,J,ks) in enumerate(cases):
    a={i:(S.Integer(1) if i in (0,I[-1]) else S.Symbol('a'+str(i))) for i in I}
    b={j:(S.Integer(1) if j in (0,J[-1]) else S.Symbol('b'+str(j))) for j in J}
    C=S.Poly(sum(a[i]*x**i for i in I)*sum(b[j]*x**j for j in J),x)
    cs={k:C.nth(k) for k in ks}
    if index==0:
        target=a[4]*(1-b[1])+b[4]*(1-a[1])+a[2]*b[2]
        assert S.expand(cs[4]-cs[5]-target)==0
        expression=str(target)
    elif index==1:
        u,v=a[5],b[4]
        target=u*(1-v)+v*(1-u)
        assert S.expand((1-cs[9]).subs({a[2]:1-v,b[7]:1-u})-target)==0
        assert cs[15]==a[2]+b[4] and cs[18]==a[5]+b[7]
        expression=str(target)
    else:
        target=(1-a[1])*b[8]+a[6]*(1-b[3])
        assert S.expand(cs[19]-cs[9]-target)==0
        assert cs[1]==a[1]+b[1]
        expression=str(target)
    outputs.append({'I':I,'J':J,'diagonals':{k:str(v) for k,v in cs.items()},
                    'positive_gap':expression,'identity_verified':True})
print(json.dumps({'sympy':S.__version__,'cases':outputs},indent=2))
