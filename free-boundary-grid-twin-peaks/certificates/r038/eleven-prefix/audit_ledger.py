import pathlib,json,math,hashlib
p=pathlib.Path(__file__).resolve().parent;d=json.loads((p/'N11_COMPLETE_LEDGER.json').read_text());source=p.parents[1]/'r037/large-side/large_side_orbit_dp.cpp'
assert hashlib.sha256(source.read_bytes()).hexdigest()==d['source_sha256']=='9c466458077713510c267d894f032dcff12e41d931b12d4ab2c841c5314fc630'
assert d['status']=='complete_exact_enumeration_authenticated' and d['controller_finished']['all_children_waited_and_reaped'] and d['controller_finished']['status']=='dispatch-complete-pending-independent-audit'
rows={l:[] for l in [6,7,8,9,10]};ids=set()
for x in d['records']:
 a=x['checkpoint'];l=a['left_size'];assert l in rows and x['wait_returncode']==x['reap_returncode']==0
 assert a['N']==11 and a['mode']=='scan-roworbit-j7' and a['status']=='completed' and a['source_sha256']==d['source_sha256'] and a['max_abs_Rj7']==0
 assert a['shard_id'] not in ids;ids.add(a['shard_id']);rows[l].append(a)
coverage=0;graphs=0
for l,rs in rows.items():
 cur=0;weight=0
 for a in sorted(rs,key=lambda x:x['orbit_start']):
  assert a['orbit_start']==cur and a['orbit_stop_exclusive']>cur
  assert a['next_orbit_rank']==a['orbit_stop_exclusive'] and a['orbits_scanned']==a['orbit_stop_exclusive']-cur
  cur=a['orbit_stop_exclusive'];weight+=a['weighted_labeled_graphs']
 assert cur==math.comb(2**(11-l)+l-1,l) and weight==2**(l*(11-l));coverage+=cur;graphs+=weight
assert len(d['records'])==613 and coverage==2501994 and graphs==1359217664
print('PASS 613 chunks; 2501994 representatives; 1359217664 labelled matrices; all R7=0')
