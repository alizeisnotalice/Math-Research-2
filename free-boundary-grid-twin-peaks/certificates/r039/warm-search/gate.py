import pathlib,subprocess,json,hashlib,os,sys,argparse
b=pathlib.Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--old',type=pathlib.Path,required=True);p.add_argument('--new',type=pathlib.Path,required=True);z=p.parse_args()
out=z.output;out.mkdir(exist_ok=False);records=[]
def run(name,cmd,good=True):
 with (out/(name+'.stdout')).open('w') as so,(out/(name+'.stderr')).open('w') as se:
  p=subprocess.Popen(cmd,stdout=so,stderr=se);rc=p.wait()
 try:os.kill(p.pid,0);absent=False
 except ProcessLookupError:absent=True
 assert absent and ((rc==0) if good else (rc!=0)),(name,rc)
 records.append(dict(name=name,pid=p.pid,returncode=rc,expected_success=good,waited_reaped=True,absent_after=absent));(out/'receipt.json').write_text(json.dumps(records,indent=2)+'\n')
old=z.old;new=z.new;verifier=b.parents[1]/'r036/finite-candidates/verify_psi.py'
for n in [3,8,32,64]:
 root=(n//2)*n+n//2;w=[root];seen={root}
 for v in w:
  for X,Y in [(v%n-1,v//n),(v%n+1,v//n),(v%n,v//n-1),(v%n,v//n+1)]:
   u=Y*n+X
   if 0<=X<n and 0<=Y<n and u not in seen:seen.add(u);w.append(u)
 inp=out/('n%d-start.word'%n);inp.write_text('n %d root %d\nword '%(n,root)+' '.join(map(str,w))+'\n')
 a=out/('n%d-old.word'%n);c=out/('n%d-warm.word'%n)
 for tag,bin,init,word in [('old',old,'bfs',a),('warm',new,'word:'+str(inp),c)]:run('n%d-%s'%(n,tag),[str(bin),str(n),str(root),'2026111739','10000',str(word),'1',init])
 aa=a.read_text().splitlines();cc=c.read_text().splitlines();assert aa[1]==cc[1] and aa[0].split(' init ')[0]==cc[0].split(' init ')[0]
 run('n%d-verify'%n,[sys.executable,'-B',str(verifier),str(c)])
# Malformed headers/permutations/roots/trailing input/one-peak all fail closed.
valid='n 3 root 4\nword 4 3 5 1 7 0 6 2 8\n'
invalid=[valid.replace('n 3','n 4'),valid.replace('root 4','root 0'),valid.replace('3 5','3 3'),valid.replace('3 5','3 9'),valid+'extra\n','n 3 root 4\nword 4 0 1 2 3 5 6 7 8\n']
for i,text in enumerate(invalid):
 inp=out/('bad%d.word'%i);inp.write_text(text);run('bad%d'%i,[str(new),'3','4','1','10',str(out/('bad%d-out.word'%i)),'1','word:'+str(inp)],False)
(out/'FROZEN_HASHES.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [b/'psi_search_warm.cpp',new,pathlib.Path(__file__)]},indent=2)+'\n');print('PASS 4 same-state trajectories and 6 malformed words',len(records),'processes')
