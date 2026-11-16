import itertools,collections,json,time,pathlib,math
start=time.time();n=3;N=n*n
adj=[sum(1<<w for w in range(N) if abs(v%n-w%n)+abs(v//n-w//n)==1) for v in range(N)]
full=collections.Counter();lift=collections.Counter();ct=collections.Counter();sumL2=0
for word in itertools.permutations(range(N)):
 seen=0;peaks=[]
 for rank,v in enumerate(word,1):
  if not adj[v]&seen:peaks.append((rank,v))
  seen|=1<<v
 if len(peaks)!=2:continue
 (unused,x),(tau,y)=peaks;d=abs(x%n-y%n)+abs(x//n-y//n)
 de=(x%n-y%n)**2+(x//n-y//n)**2
 full[(d,tau)]+=1
 ct['A2']+=1;ct['D']+=de
 if tau!=2:continue
 j=next(i for i,v in enumerate(word,1) if adj[y]>>v&1)
 L=j-2;ct['top2']+=1;ct['L']+=L;ct['LD']+=L*de;sumL2+=L*L
 for t in range(2,L+2):lift[(d,t)]+=1
assert full==lift
assert [ct[x] for x in ['A2','D','top2','L','LD']]==[143808,633600,81984,143808,633600]
assert sumL2==358080
out={'independent_direct_position_mask_checker':True,'permutations':math.factorial(N),'counts':dict(ct),'sumL2':sumL2,'all_distance_birth_cells_equal':True,'wall_seconds':time.time()-start}
pathlib.Path(__file__).with_name('result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
