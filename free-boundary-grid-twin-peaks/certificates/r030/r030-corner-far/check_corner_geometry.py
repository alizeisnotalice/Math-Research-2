from math import comb, prod
def fall(x,k): return prod(x-i for i in range(k))
for n in range(5,33):
 N=n*n; corners=(0,n-1,N-n,N-1)
 assert min(sum(abs((v%n)-(p%n))+abs((v//n)-(p//n))>=n-1 for p in corners) for v in range(N))>=2
 insert=2*comb(N-1,3)
 lhs=4*insert//2; upper=(N-3)*(N-1)*fall(N-2,2*n-3)
 assert lhs*3*fall(N,8)*fall(N-3,2*n-4)==2*upper*fall(N,8)
 assert fall(N,8)*fall(N-3,2*n-4)<=N**(2*n+4)
print({'passed':True,'scope':'geometry/count sanity n=5..32'})
