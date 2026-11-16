#!/usr/bin/env python3
"""Exact ordered-prefix MH audit and fixed-L multinomial resample-move SMC."""
from __future__ import annotations
import argparse, itertools, json, math, random
from collections import Counter
from pathlib import Path

def grid(n):
    N=n*n; adj=[0]*N
    for y in range(n):
      for x in range(n):
        k=y*n+x
        for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
          xx,yy=x+dx,y+dy
          if 0<=xx<n and 0<=yy<n: adj[k]|=1<<(yy*n+xx)
    return adj

def births(word,adj):
    seen=0; out=[]
    for t,v in enumerate(word,1):
      if not (adj[v]&seen): out.append(t)
      seen|=1<<v
    return tuple(out)

def valid(word,tau,adj):
    return births(word,adj)==((1,) if len(word)<tau else (1,tau))

def proposal(word,i,v):
    z=list(word)
    try:
      j=z.index(v); z[i],z[j]=z[j],z[i]
    except ValueError: z[i]=v
    return tuple(z)

def exact_kernel_audit(n,max_r=None,ranks=None):
    """Integer-count detailed balance and uniform-invariance check."""
    N=n*n; adj=grid(n); max_r=min(N,max_r or N); rows=[]
    for tau in range(2,N+1):
      for r in (ranks if ranks is not None else range(1,max_r+1)):
        states=[w for w in itertools.permutations(range(N),r) if valid(w,tau,adj)]
        S=set(states); outs={}; good=True
        for x in states:
          c=Counter()
          for i in range(r):
            for v in range(N):
              y=proposal(x,i,v); c[y if y in S else x]+=1
          good &= (sum(c.values())==r*N); outs[x]=c
        db=good and all(a==outs[y][x] for x,c in outs.items() for y,a in c.items())
        inv=good and all(sum(c[y] for c in outs.values())==r*N for y in states)
        rows.append({'n':n,'tau':tau,'r':r,'states':len(states),
                     'detailed_balance':bool(db),'uniform_invariant':bool(inv)})
    return rows

def allowed(word,tau,adj,N):
    rank=len(word)+1; seen=sum(1<<v for v in word)
    if rank==1: return list(range(N))
    bd=0
    for v in word: bd|=adj[v]
    bd &= ~seen
    mask=(((1<<N)-1)^seen^bd) if rank==tau else bd
    return [v for v in range(N) if mask>>v&1]

def mh_move(word,tau,L,adj,rng,N):
    for _ in range(L):
      y=proposal(word,rng.randrange(len(word)),rng.randrange(N))
      if valid(y,tau,adj): word=y
    return word

def multinomial(words,weights,rng):
    total=sum(weights)
    if not total: return None
    cut=[]; s=0.0
    for w in weights: s+=w; cut.append(s)
    out=[]
    for _ in words:
      u=rng.random()*total
      for i,z in enumerate(cut):
        if z>u: out.append(words[i]); break
    return out

def e2(word,n,adj):
    b=births(word,adj); assert len(b)==2
    x,y=word[b[0]-1],word[b[1]-1]
    return (x%n-y%n)**2+(x//n-y//n)**2

def terminal_audit(n):
    """Exact N-1-prefix check of the terminal forced set, all tau at once."""
    N=n*n; adj=grid(n); out=[]
    for tau in range(2,N+1):
      count=0; sizes=set()
      for w in itertools.permutations(range(N),N-1):
        if valid(w,tau,adj):
          count+=1; sizes.add(len(allowed(w,tau,adj,N)))
      expected={0} if tau==N else {1}
      out.append({'n':n,'tau':tau,'valid_Nminus1_prefixes':count,
                  'terminal_allowed_sizes':sorted(sizes),'expected_sizes':sorted(expected),
                  'vacuous_no_surviving_prefix':count==0,
                  'passed':(count==0 or sizes==expected)})
    return out

def smc_one(n,tau,P,L,seed):
    """GPU-matched fully-adapted schedule; count-normalizer convention."""
    N=n*n; adj=grid(n); rng=random.Random(seed); words=[tuple() for _ in range(P)]
    Z=1.0; dead=0
    for r in range(1,N+1):
      weights=[float(len(allowed(x,tau,adj,N))) for x in words]
      dead+=sum(w==0 for w in weights)
      meanw=sum(weights)/P; prior_Z=Z; Z*=meanw
      if not meanw: return 0.0,0.0,dead
      if r==N:
        # No terminal resampling. Mutate each live parent and retain final a-weight.
        sf=0.0
        for x,w in zip(words,weights):
          if w:
            a=allowed(x,tau,adj,N)
            y=x+(a[rng.randrange(len(a))],)
            sf+=w*e2(y,n,adj)
        return Z, prior_Z*sf/P, dead
      parents=multinomial(words,weights,rng)
      words=[]
      for x in parents:
        a=allowed(x,tau,adj,N)
        assert a
        y=x+(a[rng.randrange(len(a))],)
        words.append(mh_move(y,tau,L,adj,rng,N))
    raise AssertionError('unreachable')

def smc_calibration(n=3,P=1024,reps=32,Ls=(0,4,16),seed0=0xA0221800):
    ans=[]; N=n*n
    for L in Ls:
      vals=[]
      for rep in range(reps):
        z=d=0.0; dead=0
        for tau in range(2,N+1):
          a,b,c=smc_one(n,tau,P,L,seed0+1000000*L+10000*rep+tau)
          z+=a; d+=b; dead+=c
        # d/z estimates the unnormalised squared Euclidean distance; Q divides by n^2.
        vals.append((z,d/(z*n*n) if z else float('nan'),dead))
      def mean(j): return sum(x[j] for x in vals)/reps
      def se(j):
        m=mean(j)
        return math.sqrt(sum((x[j]-m)**2 for x in vals)/(reps-1)/reps)
      ans.append({'L':L,'particles':P,'replicates':reps,'seed_base':seed0,
                  'mean_A2':mean(0),'SE_A2_descriptive':se(0),
                  'mean_Q':mean(1),'SE_Q_descriptive':se(1),
                  'mean_early_dead':mean(2)})
    return ans

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out',default='resample-move-audit.json')
    ap.add_argument('--particles',type=int,default=1024); ap.add_argument('--reps',type=int,default=32)
    a=ap.parse_args()
    # n=2 all ranks; n=3 exact small-prefix state spaces, all tau and r<=4.
    # The terminal tau=N empty-set fact is handled analytically in the lemma.
    rows=exact_kernel_audit(2)+exact_kernel_audit(3,4)
    terminal=terminal_audit(2)+terminal_audit(3)
    cal=smc_calibration(P=a.particles,reps=a.reps)
    out={'method':'ordered-prefix exact detailed-balance audit; fully-adapted fixed-L SMC',
         'exact_kernel_rows':rows,
         'all_exact_checks_pass':all(x['detailed_balance'] and x['uniform_invariant'] for x in rows),
         'terminal_checks':terminal,'all_terminal_checks_pass':all(x['passed'] for x in terminal),
         'calibration':cal,'n3_exact':{'A2':143808,'Q':1100/2247},
         'caveat':'SMC values are finite numerical diagnostics, not conditional-law certification.'}
    Path(a.out).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
