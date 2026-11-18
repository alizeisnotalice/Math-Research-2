"""Join independently replayed A-candidate bins after checking the entire reference grid and both conditional signs."""
import argparse,json,sys
from hashlib import sha256
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx


def run(base,files,out):
    ctx.prec=128;chosen={};hashes={};baseline=None;domain=json.loads((base/'sigma-154-a-grade-reference-domain.json').read_bytes());assert arb(domain['min_signed_gap_to_either_end'])>0
    for name in files:
        path=base/name;raw=path.read_bytes();q=json.loads(raw);settings=q['parameters'];assert settings['line']==domain['line']==['.429','.009'] and settings['sigma_degrees']==domain['sigma']=='1.54' and settings['physical_parameters']==['.273','.306','.0102'];hashes[name]=sha256(raw).hexdigest()
        assert all(sha256((base/k).read_bytes()).hexdigest()==v for k,v in q['source_hashes'].items())
        if baseline is None:baseline=q['source_hashes']
        assert q['source_hashes']==baseline and q['bits']>=128
        assert len(q['rows'])==len(settings['bins'])
        for row,(j,mode,params) in zip(q['rows'],settings['bins'],strict=True):
            assert row['bin']==j and 0<=j<70
            if not row['positive']:continue
            ranges=row['reachable'];assert ranges==([[0,6]] if j==69 else [[0,75]])
            classes=row['classes'];assert len(classes)==(2 if mode=='split' else 1)
            if mode=='split':
                l,a,b,r,tau=map(F,params)
                assert [(list(map(F,s['bracket'])),s['anchor_comparison_sign'],F(s['anchor_degrees']),F(s['anchor_threshold'])) for s in classes]==[([l,a],-1,a,-tau),([b,r],1,a,-tau)]
            else:
                assert mode=='direct' and list(map(F,classes[0]['bracket']))==list(map(F,params[:2])) and classes[0]['pose']==params[2:]
            assert all(s['reference_bin']==j and s['subcells']==ranges[0] for s in classes)
            root=min(arb(s['root_margin_lower']).lower() for s in classes);clear=min(arb(s['placement_margin_lower']).lower() for s in classes);assert min(root,clear)>0
            if j not in chosen or clear>chosen[j][0]:chosen[j]=(clear,root,name,row)
    assert set(chosen)==set(range(70)),sorted(set(range(70))-set(chosen));low=min((x[0],j) for j,x in chosen.items());root=min((x[1],j) for j,x in chosen.items());count=sum(b-a for x in chosen.values() for a,b in x[3]['reachable']);assert count==5181
    doc={'sigma_degrees':'1.54','AH_line':['.429','.009'],'EH_relocation_row_depth':'.000555','physical_gate':domain,'certificates_sha256':hashes,'frozen_sources_sha256':baseline,'aggregate_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'number_reference_bins':len(chosen),'reachable_reference_subcells':count,'root_margin_lower':str(root[0]),'root_margin_bin':root[1],'placement_margin_lower':str(low[0]),'placement_margin_bin':low[1],'placement_rows_above_conservative_EH_transfer':bool(low[0]>arb('.00000009488')),'number_placed_leaves':sum(s['placed_leaves'] for x in chosen.values() for s in x[3]['classes']),'number_excluded_leaves':sum(s['root_free_leaves'] for x in chosen.values() for s in x[3]['classes']),'chosen_bins':[{ 'bin':j,'certificate':chosen[j][2],'root_margin':str(chosen[j][1]),'placement_margin':str(chosen[j][0])} for j in range(70)],'python':sys.version};out.write_text(json.dumps(doc,indent=2)+'\n');print(root,low,count)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('output',type=Path);ap.add_argument('certificates',nargs='+');args=ap.parse_args();run(args.archive.resolve(),args.certificates,args.output)
