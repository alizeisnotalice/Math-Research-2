"""Eliminate an isolated linear variable even from a nonlinear equation."""
from modular_core_29_30 import substitution


def forward_eliminate(R,variables,equations,track):
    pairs=[]
    n=len(equations)
    for k,g in enumerate(equations):
        if g:
            history=[R.one if i==k else R.zero for i in range(n)] if track else None
            pairs.append((g,history))
    while True:
        unit=next((pair for pair in pairs if pair[0].is_ground),None)
        if unit is not None:
            return [unit]
        best=None
        for k,(g,_) in enumerate(pairs):
            for i,x in enumerate(variables[:-1]):
                pivot=x.LM
                if pivot not in g or any(m[i] and m!=pivot for m in g):
                    continue
                degrees=[sum(m) for m in g if m!=pivot]
                degree=max(degrees,default=0)
                trial_degree=max((sum(m)+m[i]*(degree-1) for r,(f,_) in enumerate(pairs) if r!=k for m in f),default=0)
                if degree>1 and trial_degree>8:
                    continue
                rank=(degree>1,trial_degree,degree,len(g),-i,k)
                if best is None or rank<best[0]:
                    best=(rank,k,i,g[pivot])
        if best is None:
            return pairs
        _,k,i,coef=best
        g,history=pairs.pop(k)
        inv=R.domain.one/coef
        x=variables[i]
        value=x-g*inv
        assert all(m[i]==0 for m in value)
        pivot_history=[f*inv for f in history] if track else None
        revised=[]
        for f,source in pairs:
            if all(m[i]==0 for m in f):
                revised.append((f,source))
                continue
            reduced,quotient=substitution(R,f,i,value,x)
            if reduced:
                new_source=[a-quotient*b for a,b in zip(source,pivot_history,strict=True)] if track else None
                revised.append((reduced,new_source))
        pairs=revised
