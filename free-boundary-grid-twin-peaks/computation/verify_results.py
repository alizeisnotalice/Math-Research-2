"""Independent consistency audit of archived finite certificates and witnesses."""
from fractions import Fraction
import json
from pathlib import Path
from search import independent_peaks

p=Path(__file__).parent/'results'
load=lambda name:json.loads((p/name).read_text())
a={r['n']:r for r in load('exact-r001.json')['results']}
pairs=load('pairs-r001.json')['results']
m=load('moment-n5-r001.json')
for r in pairs:
 n=r['n'];assert sum(r['distance_counts'].values())==r['A2']
 if n==4:assert r['distance_counts']==a[n]['distance_counts'] and r['squared_euclidean_sum']==a[n]['squared_euclidean_sum']
 if n==5:assert r['A2']==m['A2'] and r['squared_euclidean_sum']==m['squared_euclidean_sum']
for r in load('cuts-r001.json')['results']:
 counts=r['intersections'];union=sum((1 if len(k.split(','))%2 else -1)*v for k,v in counts.items())
 assert union==r['union_count'] and Fraction(r['union_ratio'])==Fraction(union,r['A2'])
 assert 0<Fraction(r['second_moment_bound'])<=Fraction(r['union_ratio'])<=1
for name in ('walks-r001.json','walks-reflect-r001.json','walks-block-r001.json'):
 for r in load(name)['results']:
  n=r['n'];labels=r['best_labels'];assert sorted(labels)==list(range(1,n*n+1))
  pp=independent_peaks(labels,n);assert len(pp)==2
  u,v=pp;assert abs(u%n-v%n)+abs(u//n-v//n)==r['best_distance']
  assert sum(r['distance_histogram'])==r['samples']
  assert all(x+y==z for x,y,z in zip(r['first_half_histogram'],r['second_half_histogram'],r['distance_histogram']))
for r in load('search-r001.json')['runs']:
 if r['best']:
  n=r['n'];pp=independent_peaks(r['best']['labels_row_major'],n);assert len(pp)==2
  u,v=pp;assert abs(u%n-v%n)+abs(u//n-v//n)==r['best']['distance']
tails=[json.loads(line) for line in (p/'n5-tail.jsonl').read_text().splitlines()]
assert int(tails[0]['tail_count'])==m['A2']
assert all(int(x['tail_count'])>=int(y['tail_count']) for x,y in zip(tails,tails[1:]))
assert load('insertion-r001.json')['pairwise_match']
g={r['n']:r for r in load('geodesic-r001.json')['results']}
g4=load('geodesic-cpp-n4.json')
assert [(r['d'],r['numerator']) for r in g[4]['rows']]==[(r['d'],r['numerator']) for r in g4['rows']]
for n,r in g.items():
 assert len(r['rows'])==2*n-1
 for row in r['rows']:
  if row['d']>=2:assert Fraction(row['p_d'])<=Fraction(1,2**(row['d']-2))
print(json.dumps({'status':'passed','checks':['independent counts and squared moments','all cut intersections and inclusion-exclusion','all saved construction witnesses','monotone merge tails','insertion pairwise enumeration','full geodesic prefixes cross-language']}))
