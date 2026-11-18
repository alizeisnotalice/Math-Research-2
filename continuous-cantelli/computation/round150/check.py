import json,mpmath as m
m.mp.dps=48
p=lambda s,x: m.exp(-x*x/(2*s))/m.sqrt(2*m.pi*s)
v,w,M=map(m.mpf,['1.1','.4','3.2'])
rows=[]
for t in map(m.mpf,['.1','1','6.5']):
    def kernel(n):
        u=v*M-w*n
        if n<=0 or u<=0:return 0
        return (v-w)*M*p(u,n-M)/u*p(w*n,t-n)
    mass=m.quad(kernel,[0,M/2,M,2*M,v*M/w])
    diff=abs(mass-p(v*M,t-M))
    assert diff<m.mpf('1e-38')
    rows.append([str(t),m.nstr(diff,10)])
print(json.dumps(dict(mpmath=m.__version__,digits=m.mp.dps,upper_distance_reciprocal=str(v),lower_distance_reciprocal=str(w),m=str(M),checks=rows)))
