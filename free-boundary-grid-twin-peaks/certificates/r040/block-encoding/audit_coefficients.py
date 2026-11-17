from math import comb,factorial
count=0
for K in range(2,101):
 for b in range(1,K):
  assert sum(comb(q-1,b) for q in range(b+1,K+1))==comb(K,b+1)
  count+=1
 assert sum(2**(b-1)*comb(K-1,b) for b in range(1,K))==(3**(K-1)-1)//2
 assert sum(2**(b-1)*comb(K,b+1) for b in range(1,K))==(3**K-2*K-1)//4
 assert (3**(K-1)-1)//2<=3**(K-1)-2**(K-1)
for N in range(3,101):
 for k in range(1,N-1):
  total=sum(2**(b-1)*comb(N-1,b+1) for b in range(1,k+1))
  assert total*factorial(k+1)<=2**k*N**(k+1)
print('PASS',count,'hockey-stick coefficients; all q sums, improvements, and small-branch majorants for stated ranges')
