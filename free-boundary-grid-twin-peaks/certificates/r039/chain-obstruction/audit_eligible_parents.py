import json
out=[]
for n in [4,5,8,16,32,64,128]:
 N=n*n;nb=lambda v:[Y*n+X for X,Y in [(v%n-1,v//n),(v%n+1,v//n),(v%n,v//n-1),(v%n,v//n+1)] if 0<=X<n and 0<=Y<n]
 early=[[u for u in nb(v) if u<v] for v in range(N)];assert all(early[v] for v in range(1,N))
 parents={u[0] for u in early if len(u)==1};eligible=parents-({0}|set(nb(0)));expected={x for x in range(2,n-1)}|{y*n for y in range(2,n-1)}
 assert eligible==expected and len(eligible)==2*(n-3)
 out.append({'n':n,'eligible_sole_child_parents':len(eligible),'universal_chain_score_upper_bound':8*(n-3)**3,'N_squared':N*N,'exact_set_match':True})
print(json.dumps(out,indent=2))
