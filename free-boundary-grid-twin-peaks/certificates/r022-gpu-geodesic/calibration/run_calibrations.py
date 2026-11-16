import json,subprocess,pathlib
cfg=json.load(open('launch.json'))
for r in cfg['runs']:
 m=r['metadata'];name=r['name'];raw=name+'.jsonl';assert not pathlib.Path(raw).exists()
 subprocess.run(['./geodesic_prefix_r022',str(m['n']),str(m['particles']),str(m['replicates']),str(m['tau_lo']),str(m['tau_hi']),raw,m['domain'],'strat','16'],check=True)
