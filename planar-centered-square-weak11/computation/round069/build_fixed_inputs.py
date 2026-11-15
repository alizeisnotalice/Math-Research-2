#!/usr/bin/env python3
"""Extract exactly six pinned maximum-block states plus two structural controls."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path

T='planar-centered-square-weak11'
PINS={
 'computation/frontier063/structural-inputs.json':'20ffbb4a3e444d1b72b8fa6e3bfcafaf839f71144bd6298a9419b28af786ce41',
 'computation/results/round068-max-tree/classification.json':'2d4ef1623ed67aa661a09f606b8222d4ee8e9d3d559f9a7c13b7b35898203f53',
 'computation/inputs/r068-reachable-surplus-six-source.json':'98c4fdc86206906a2ccf7302c56746e340e0d7141d3eab288ce449170730a307',
 'computation/inputs/r068-support-tree-payment.json':'25ccd96bb2e35a37c78ce680a12654a67d25db8adcb87a270d4c22b47cce7aa4'}
STATES=[('s00-c00-n16',48623,'91/8'),('s00-c01-n17',114159,'197/16'),
        ('s10-c00-n16',24575,'77/8'),('s10-c01-n17',90111,'189/16'),
        ('s10-c02-n18',90111,'55/4'),('s10-c02-n18',179711,'55/4')]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-archive-root',type=Path,required=True)
    ap.add_argument('--out-root',type=Path,required=True);args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    data={}
    for rel,h in PINS.items():
        p=args.source_archive_root/T/rel;assert sha(p)==h,(rel,sha(p));data[rel]=json.loads(p.read_text())
    old={r['id']:r for s in data['computation/frontier063/structural-inputs.json']['strata'] for r in s['inputs']}
    cls=data['computation/results/round068-max-tree/classification.json']['blocks'];jobs=[]
    remaining={(r['input_id'],r['mask'],r['L']) for r in cls if r['priority_class']=='remaining'}
    assert remaining==set(STATES)
    def save(name,payload):
        rel=f'{T}/computation/inputs/r069/{name}.json';p=args.out_root/rel;p.parent.mkdir(parents=True,exist_ok=True)
        payload['id']=name;payload['builder_sha256']=sha(Path(__file__))
        p.write_text(json.dumps(payload,indent=2)+'\n')
        jobs.append({'job_id':name,'input':rel,'input_sha256':sha(p),
                     'output':f'{T}/computation/results/round069-state/{name}.json'})
    for identity,mask,L in STATES:
        r=old[identity];row=next(c for c in cls if c['input_id']==identity and c['mask']==mask)
        assert F(row['L'])==F(L)
        maxima=sorted(c['mask'] for c in cls if c['input_id']==identity)
        payload={'source_input_id':identity,'points':[[str(F(x,4)),str(F(y,4))] for x,y in r['points_quarters']],
                 'initial_masses':[str(F(v,16)) for v in r['mass_sixteenths']],
                 'S_mask':mask,'expected_L':L,'expected_maximum_masks':maxima,
                 'provenance':{'input_archive':f'{T}/computation/frontier063/structural-inputs.json',
                    'input_sha256':PINS['computation/frontier063/structural-inputs.json'],
                    'classification_archive':f'{T}/computation/results/round068-max-tree/classification.json',
                    'classification_sha256':PINS['computation/results/round068-max-tree/classification.json']},
                 'scope':'One previously unclassified global maximum state; original rational coordinates/masses retained.'}
        assert sum(F(v) for i,v in enumerate(payload['initial_masses']) if mask>>i&1)==F(L)
        save(f'{identity}-S{mask}',payload)
    six=data['computation/inputs/r068-reachable-surplus-six-source.json']
    assert six['points'][:4]==[['3/10','3/10'],['6/5','6/5'],['3/10','6/5'],['6/5','3/10']]
    save('r068-final-four-corner-six-source-S15',{'source_input_id':six['id'],
        'points':six['points'],'initial_masses':six['initial_masses'],'S_mask':15,'expected_L':'1','expected_maximum_masks':[15],
        'provenance':{'archive':f'{T}/computation/inputs/r068-reachable-surplus-six-source.json',
            'sha256':PINS['computation/inputs/r068-reachable-surplus-six-source.json']},
        'scope':'Final four-corner six-source control, not the discarded centered variant.'})
    tree=next(r for r in data['computation/inputs/r068-support-tree-payment.json']['examples'] if r['id']=='two-row-tree-method-barrier')
    assert len(tree['points'])==9 and tree['S_count']==8
    save('r068-tree-barrier-nine-source-S255',{'source_input_id':tree['id'],'points':tree['points'],
        'initial_masses':tree['masses'],'S_mask':255,'expected_L':'1','expected_maximum_masks':[255],
        'provenance':{'archive':f'{T}/computation/inputs/r068-support-tree-payment.json',
            'sha256':PINS['computation/inputs/r068-support-tree-payment.json'],'selected_example':tree['id']},
        'scope':'Eight-point maximum S plus one outside source; tree criterion barrier control.'})
    spec={'schema':'r069-fixed-eight-state-spec-v1','jobs':jobs,'bits':64,
          'geometry_source':f'{T}/computation/round064/mass_polytope_variation.py',
          'geometry_sha256':'4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246',
          'rectclip_source':f'{T}/computation/round066/erosion_ring_diagnostic.py',
          'rectclip_sha256':'4dc4fb494736c40334359f9950f052ee2e4a53b503883ffbc3f9bc0da69b1361',
          'K':'[bboxS.xmax-sqrt(3L)/2,bboxS.xmin+sqrt(3L)/2] x corresponding y interval',
          'source_pins':PINS,'scope':'Exactly six historical maximum states in five inputs plus final six-source and tree-nine-source controls; no random/new family.'}
    sp=args.out_root/T/'computation/round069/fixed-eight-job-spec.json';sp.parent.mkdir(parents=True,exist_ok=True)
    sp.write_text(json.dumps(spec,indent=2)+'\n');print(json.dumps({'jobs':len(jobs),'spec':str(sp),'sha256':sha(sp)}))


if __name__=='__main__':main()
