#!/usr/bin/env python3
"""Final R66 review provenance/receipt binding; no producer recomputation."""
import hashlib,json,sys
from pathlib import Path
from fractions import Fraction as F
HERE=Path(__file__).resolve().parent;TOP=HERE.parent
ROOT=Path('/Volumes/data/CodexProjects/T2');FR=ROOT/'.local/compute/frontier066'
CP=ROOT/'.local/workspaces/pcw-position-slab-impl/planar-centered-square-weak11'
AP=ROOT/'.local/workspaces/pcw-astra-covering/planar-centered-square-weak11'
ER=CP/'computation/results/round066-erosion'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def pin(p):return dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size)
BS='211d33e21e47e6e14290e4ddbc3b3a36e247131de1bc53114091b7766d0fe8b6'
ES='4dc4fb494736c40334359f9950f052ee2e4a53b503883ffbc3f9bc0da69b1361'
MS='4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246'
PROOFS={
'r066-layered-high-competition-budget.md':'6c1530bf49240e446d6c8f4e626ddf5d9e23566bbd7fe8fd34d2b3af3de5499a',
'r066-erosion-rim-and-repeated-capacity.md':'91761feabcd30422ba52f03f9b35ab5f02eefaa082dded717685391fb7c80cb5',
'r066-five-source-shape-loss-obstruction.md':'7ea8e042c6d501a4a2fa49e511411bacc6118c99e6a598855bab442674743a2a',
'r066-three-source-fat-envelope-obstruction.md':'6c17a3a9c98a92a54b240f962d6db25aa3c706448bacbca67920ec00733f805b'}
def main():
    assert Path('/Volumes/data').is_mount();sources=set();checks=[]
    for n,h in PROOFS.items():p=AP/'proofs'/n;assert sha(p)==h;sources.add(p)
    for p,h in [(FR/'theta_competitor_budget.py',BS),(CP/'computation/round066/theta_competitor_budget.py',BS),(CP/'computation/round066/erosion_ring_diagnostic.py',ES),(ROOT/'.local/compute/r064-lp-pilot/mass_polytope_variation.py',MS)]:assert sha(p)==h;sources.add(p)
    # Bind all independent-result provenance pins and own source hashes.
    for name in ('round066-independent-grid-replay.json','round066-independent-critical-controls.json'):
        out=read(HERE/'results'/name)
        for v in out['source_pins']:
            p=Path(v['path']);assert sha(p)==v['sha256']and p.stat().st_size==v['bytes'];sources.add(p)
    assert read(HERE/'results/round066-independent-grid-replay.json')['script_sha256']==sha(HERE/'review_r066_grid_budgets.py')
    assert read(HERE/'results/round066-independent-critical-controls.json')['script_sha256']==sha(HERE/'review_r066_critical_controls.py')
    def receipt(p):
        d=read(p);assert d['exit_status']==0 and d['owned_pid_absent']is True;sources.add(p);return d
    specs=[(FR/'r64-final-theta.json',ER/'r64-final-erosion.json',FR/'r64-pilot-receipt.json',ER/'r64-final.receipt.json')]
    specs += [(FR/f'stratum-{i:02}.json',ER/f'stratum-{i:02}-erosion.json',FR/f'stratum-{i:02}.receipt.json',ER/f'stratum-{i:02}.receipt.json')for i in range(12)]
    for bp,ep,brp,erp in specs:
        b,e=read(bp),read(ep);br,er=receipt(brp),receipt(erp)
        assert b['source_sha256']==BS and b['geometry_module_sha256']==MS
        inp=Path(b['input_path']);assert sha(inp)==b['input_sha256'];sources.add(inp)
        assert e['source_sha256']==ES and e['geometry_sha256']==MS and e['budget_result_sha256']==sha(bp)
        assert er['source_sha256']==ES and er['geometry_sha256']==MS and er['budget_result_sha256']==sha(bp)and er['output_sha256']==sha(ep)
        assert br['command'][br['command'].index('--out')+1]==str(bp)
        checks.append(dict(budget=str(bp),erosion=str(ep),cases=len(b['cases']),blocks=b['evaluated_proper_blocks'],source_and_input_chain=True,receipts_saved_exit0_absent=True))
    batch=read(FR/'batch-receipt.json');assert batch['status']=='PASS'and(batch['total_cases'],batch['total_blocks'])==(72,199);sources.add(FR/'batch-receipt.json')
    for p in (ER/'batch-pilot.receipt.json',ER/'batch-strata.receipt.json'):assert read(p)['status']=='PASS';sources.add(p)
    idx=read(CP/'computation/round066/actual-erosion-file-index.json');copies=0
    for v in idx['files']:
        a,b=Path(v['original']),Path(v['copy']);assert sha(a)==sha(b)==v['sha256']and b.stat().st_size==v['bytes'];copies+=1
    sources.add(CP/'computation/round066/actual-erosion-file-index.json')
    for role in ('five-source-control','three-source-fat-control'):
        p=ER/role;r=read(p/'receipt.json');b=read(p/'budget.json');e=read(p/'erosion.json')
        assert r['budget_output_sha256']==sha(p/'budget.json')and r['erosion_output_sha256']==sha(p/'erosion.json')
        assert e['budget_result_sha256']==sha(p/'budget.json')
        assert r['module_sha256']==MS
        for c in r['receipts']:assert c['exit_status']==0 and c['owned_pid_absent']is True
        if role=='three-source-fat-control':
            assert r['input_sha256']==b['input_sha256']==sha(p/'input.json')
            assert r['budget_source_sha256']==BS and r['erosion_source_sha256']==ES
            a=read(AP/'computation/inputs/r066-three-source-fat-envelope.json');adapt=read(p/'input.json')
            assert a['points']==adapt['points']and list(map(F,a['masses']))==list(map(F,adapt['initial_masses']))and a['S_mask']==adapt['S_mask']==3
            original=Path('/Volumes/data/CodexTemp/T2/round066/three-source-fat-control-20261115')
            for n in ('input.json','budget.json','erosion.json','receipt.json'):assert sha(original/n)==sha(p/n)
        else:assert r['budget_sha256']==BS and r['erosion_sha256']==ES
        sources.add(p/'receipt.json')
    # External actual executions: validate records, do not claim independent rerun.
    w=receipt(FR/'wolfram-v2-receipt.json');wp=FR/'r066-five-source-independent-v2.wl'
    assert w['source_sha256']==sha(wp)=='1085a949bf4acdd0abef07857ec33e3aa4196bdcf79f5e2aec051cce89b51da0'
    wd=read(FR/'r066-five-source-wolfram-v2.json');assert wd['status']=='PASS_EXACT_ALGEBRAIC_CHECKS'and all(wd['checks'].values())
    sources.update((wp,FR/'r066-five-source-wolfram-v2.json'))
    fp=FR/'check_fat_repeated_capacity.py';fd=read(FR/'check_fat_repeated_capacity.json')
    assert fd['script_sha256']==sha(fp)=='524d2884d955d85f71237e5594610ee95beef4859e819f093650117120630b91'
    assert fd['status']=='PASS_EXACT_FINITE_FAMILY_CHECKS'and fd['total_closed_bboxes']==1009644 and all(fd['owned_workers_absent'].values())
    sources.update((fp,FR/'check_fat_repeated_capacity.json'))
    for n in ('computation/inputs/r066-three-source-fat-envelope.json','computation/results/round066/three-source-fat-arithmetic.json'):sources.add(AP/n)
    result=dict(status='PASS_R66_FINAL_BINDING',python=sys.version,proof_pins=PROOFS,
        original_budget_erosion_pairs=checks,bytewise_erosion_inventory_copies=copies,
        separate_controls=['five-source-thin','three-source-fat'],source_pins=[pin(p)for p in sorted(sources)],
        scope='Saved source/receipt/hash-chain binding plus completed independent exact grid replay. External Wolfram and 12-N census not rerun. No general Delta/C2 result.')
    dest=HERE/'results/round066-independent-final-binding.json';dest.write_text(json.dumps(result,indent=2)+'\n')
    owned=[HERE/'review_r066_grid_budgets.py',HERE/'review_r066_critical_controls.py',Path(__file__),HERE/'results/round066-independent-grid-replay.json',HERE/'results/round066-independent-critical-controls.json',dest,TOP/'reviews/r066-independent-budget-and-obstruction-review.md']
    manifest=dict(round=66,status='FROZEN_INDEPENDENT_AUDIT_DELIVERY',files=[pin(p.resolve())for p in owned],final_proof_pins=PROOFS,
        scope='201-block original replay and separate new fat control; exact identities/erosion/obstruction proof review, finite arithmetic evidence, final provenance binding. No producer edits, MATLAB, commit or push.')
    mp=HERE/'results/round066-independent-manifest.json';mp.write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],source_pin_count=len(sources),artifact_count=len(owned),manifest_sha256=sha(mp),manifest_path=str(mp)),indent=2))
if __name__=='__main__':main()
