import json,subprocess,pathlib
for r in json.load(open('launch.json'))['runs']:
 m=r['metadata'];name=r['name'];raw=name+'.jsonl';assert not pathlib.Path(raw).exists()
 subprocess.run(['./geodesic_prefix_r022',str(m['n']),str(m['particles']),str(m['replicates']),str(m['tau_lo']),str(m['tau_hi']),raw,m['domain'],'strat','16'],check=True)
 subprocess.run(['python3','summarize_onepeak_geodesic.py','--input',raw,'--sidecar',name+'.sidecar.json','--source','geodesic_prefix_r022.swift','--out',name+'.summary.json'],check=True)
