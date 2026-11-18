import json,re,datetime
from decimal import Decimal
from pathlib import Path
def read(n):return json.loads(Path(n).read_text())
diff=read('offset.json');o0=max(v[0] for v in diff);o1=min(v[1] for v in diff)
assert o0<=o1
def first_common(name):
 rows=read(name); first={}
 for row in rows:
  i=row[1]
  if i not in first or row[2]<first[i][2]:first[i]=row
 return max(r[2] for r in first.values()),min(r[3] for r in first.values()),len(first)
s0,e0,n0=first_common('local.json');s1,e1,n1=first_common('mac.json')
s=max(s0,s1-o0);e=min(e0,e1-o1)
samples=[]
for time,nat in re.findall(r'utc=(\S+) batches=\d+ tested=\d+ gpu_seconds=(\S+)',Path('localgpu.txt').read_text()):
 t=datetime.datetime.fromisoformat(time.replace('Z','+00:00')).timestamp()
 if s+.002<t<e-.002:samples.append((t,float(nat)))
local=samples[-1][1]-samples[0][1]
remote=0
for a,b,nat in re.findall(r'utc0=(\S+) utc1=(\S+) gpu=(\S+)',Path('macgpu.txt').read_text()):
 if float(a)-o1>=s and float(b)-o0<=e:remote+=float(nat)
over=local+remote-(e-s)
rows=read('local.json')+read('mac.json');rows.sort(key=lambda r:Decimal(r[5]))
assert all(Decimal(r[6])==Decimal(q[5]) for r,q in zip(rows,rows[1:]))
upper=max(Decimal(r[-1][1]) for r in rows);assert upper<0 and over>0 and (n0,n1)==(12,15)
print(json.dumps(dict(workers=[n0,n1],offset=[o0,o1],assumption='stable remote-minus-local clock offset within measured bracket',local_unix_interval=[s,e],cpu_window_seconds=e-s,local_gpu_seconds_lower=local,mac_gpu_seconds_lower=remote,joint_device_command_seconds_lower=over,boxes=len(rows),left=rows[0][5],right=rows[-1][6],upper=str(upper)),indent=2))
