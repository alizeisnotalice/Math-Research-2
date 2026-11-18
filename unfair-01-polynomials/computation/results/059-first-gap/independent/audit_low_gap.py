"""Third Python/SymPy replay of the first-gap polynomial certificates."""
from argparse import ArgumentParser
from pathlib import Path
from hashlib import sha256
import json
import sys
import sympy as S


def fail_if(b):
    if b:
        raise RuntimeError('identity verification failed')


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args()
    t=S.symbols('t')
    m=1-t+t*t
    v=1-(1-t)*m
    u=1-t*v
    f=2*t**3-4*t**2+5*t-1
    fail_if(S.expand((1-t)*u+m*v-1-t*(t-1)*f)!=0)
    values={str(c):str(f.subs(t,c)) for c in (-1,1,S.Rational(-1,2),S.Rational(1,2))}
    fail_if(any(value=='0' for value in values.values()) or not S.Poly(f,t).is_irreducible)
    n=1-t*m
    p=1-(1-t)*n
    w=S.expand(1-(1-t)*p-m*n)
    fail_if(S.expand(w-t*t*m)!=0)
    H=t**6-4*t**5+6*t**4-7*t**3+4*t**2-2*t+1
    fail_if(S.expand(t*w+p*m-1+H)!=0)
    h2=S.Poly(H,t,modulus=2)
    powers=[pow(2,k,9) for k in range(1,7)]
    fail_if(powers!=[2,4,8,7,5,1] or not h2.is_irreducible)
    t_bounds=(S.Rational(7,10),S.Rational(3,4))
    H_poly=S.Poly(H,t)
    fail_if(H_poly.count_roots(0,1)!=1 or H_poly.count_roots(*t_bounds)!=1)
    A=[1,t,0,n,p,0,0,0,0]
    B=[1,1-t,m,0,0,w,0,0,0]
    c=[S.rem(S.expand(sum(A[i]*B[k-i] for i in range(k+1))),H,t) for k in range(9)]
    fail_if(c[:8]!=[1]*7+[0] or S.rem(c[8]-n*w,H,t)!=0)
    data={'status':'verified','sympy':S.__version__,'python':sys.version,
          'source_sha256':sha256(a.source.read_bytes()).hexdigest(),'c6_cubic_hypothetical_rational_roots':values,
          'C7_min_polynomial_mod2':str(h2.as_expr()),'powers_of_two_mod_nine':powers,
          'unique_real_root_t_in_unit_interval':True,
          'H_at_bounds':[str(H.subs(t,x)) for x in t_bounds],
          'local_product_coefficients_through_degree_eight':[str(x) for x in c]}
    with a.out.open('x') as f:
        f.write(json.dumps(data,sort_keys=True,indent=2)+'\n')
    print(json.dumps(data,sort_keys=True))


if __name__=='__main__':
    main()
