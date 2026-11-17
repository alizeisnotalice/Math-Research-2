"""Compare all independent screens and all AB/HGQ trees, then compact."""
import argparse,base64,gzip,hashlib,json,os,time
from pathlib import Path
from collections import Counter
def sha(b):return hashlib.sha256(b).hexdigest()
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def normalized_independent_counts(counts,codes):
    hist=Counter(codes)
    for k in range(4,8):assert counts.get(f'reject{k}',0)==hist[k]
    return Counter({k.replace('axis_','axis'):v for k,v in counts.items() if not k.startswith('reject')})
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();t=a.root;r=a.out;r.mkdir(parents=True,exist_ok=True);(r/'blocks').mkdir(exist_ok=True);start=time.monotonic()
    parent=t/'round050-parent/s2-screen-complete';ind=t/'agents/algebraic_global/round050-s2/independent-screen/full';tree=t/'round050-parent/s2-hgq-complete';abroot=t/'agents/section_eigen39/round050-s2-bfs'
    raw=(abroot/'bfs_s2_all359_depth96.json').read_bytes();ab=json.loads(raw);absha=sha(raw);inp=json.loads((abroot/'inputs_s2_all359_frozen.json').read_bytes());assert ab['input_sha256']==sha((abroot/'inputs_s2_all359_frozen.json').read_bytes())
    abmap={tuple(z['row']):z for z in ab['records']};assert len(abmap)==len(ab['records'])==9775
    assert len(inp['rows'])==9775 and set(map(tuple,inp['rows']))==set(abmap)
    assert {z['input_index'] for z in ab['records']}==set(range(9775))
    for z in ab['records']:assert z['row']==inp['rows'][z['input_index']]
    coremap={z['index']:z for z in inp['source_files']};assert set(coremap)==set(range(359))
    prefix=t/'agents/algebraic_global/round050-s2/s2_low_prefix_result.json';pcode=t/'round050-parent/s2_full_screen.py';qcode=t/'agents/algebraic_global/round050-s2/independent-screen/s2_full_mother_screen.py'
    binding={'prefix':sha(prefix.read_bytes()),'parent':sha(pcode.read_bytes()),'independent':sha(qcode.read_bytes())};states=json.loads(prefix.read_text());cores=sorted({tuple(z[k] for k in ('T','N','Delta','b','c')) for z in states['survivors']});assert len(cores)==359
    blocks=[];pending=[];seen=set();total=Counter();depths=Counter();layers=0;maxlive=0
    for i in range(359):
        qp=ind/f'core-{i:03}.json'
        while not qp.exists() and time.monotonic()-start<14400:time.sleep(2)
        if not qp.exists():raise RuntimeError('independent screen incomplete; no complete report')
        pp=parent/f'core-{i:03}.json';pb=pp.read_bytes();qb=qp.read_bytes();p=json.loads(pb);q=json.loads(qb);h=json.loads((tree/f'core-{i:03}.json').read_bytes())
        assert p['index']==q['index']==h['index']==i and p['status']==q['status']=='complete'
        assert p['core']==q['core']==list(cores[i]);assert p['binding']=={'source':binding['parent'],'prefix':binding['prefix']}
        assert q['binding']=={'source_sha256':binding['independent'],'prefix_sha256':binding['prefix']}
        assert coremap[i]['sha256']==h['screen_sha256']==sha(pb)
        packed=(parent/f'core-{i:03}-codes.bin.gz').read_bytes();codes=gzip.decompress(packed);other=(ind/f'core-{i:03}-codes.bin.gz').read_bytes()
        assert sha(packed)==p['ledger_sha256'] and sha(other)==q['gzip_code_sha256']
        assert codes==gzip.decompress(other) and len(codes)==105315 and sha(codes)==p['raw_code_sha256']==q['raw_code_sha256']
        assert Counter(p['counts'])==normalized_independent_counts(q['counts'],codes) and p['survivors']==q['survivors'];total.update(p['counts'])
        hm={tuple(z['row']):z for z in h['records']};assert len(hm)==len(h['records'])==len(p['survivors']) and set(hm)==set(map(tuple,p['survivors']))
        cert=[]
        for row in p['survivors']:
            key=tuple(row);assert key not in seen;seen.add(key);x,y=abmap[key],hm[key]
            assert x['status']==y['status']=='dead_all_prefixes' and x['depth_reached']==y['depth_reached'] and x['trace']==y['trace']
            death=y['depth_reached'];assert death==len(y['trace']) and y['trace'][-1]['live']==0
            ml=max(z['live'] for z in y['trace']);maxlive=max(maxlive,ml);layers+=death;depths[death]+=1
            cert.append({'row':row,'death':death,'max_live':ml,'trace_sha256':sha(canonical(y['trace']))})
        pending.append({'index':i,'parent_utf8':pb.decode(),'parent_sha256':sha(pb),'independent_utf8':qb.decode(),'independent_sha256':sha(qb),'codes_base64':base64.b64encode(packed).decode(),'codes_gzip_sha256':sha(packed),'records':cert})
        if len(pending)==32 or i==358:
            name=f'blocks/core-{i-len(pending)+1:03}-{i:03}.jsonl.gz';raw=b''.join(canonical(z)+b'\n' for z in pending);z=gzip.compress(raw,mtime=0);(r/name).write_bytes(z);blocks.append({'path':name,'bytes':len(z),'sha256':sha(z),'start':i-len(pending)+1,'count':len(pending)});pending=[]
    assert seen==set(abmap) and len(seen)==9775 and layers==76172 and max(depths)==67 and maxlive==4
    report={'status':'passed_all_359_screens_and_9775_exact_tree_traces','cores':359,'raw_mothers':37808085,'mothers':9775,'layers':layers,'max_death':max(depths),'max_live':maxlive,'death_histogram':dict(sorted(depths.items())),'counts':dict(total),'binding':binding,'ab_result_sha256':absha,'elapsed_seconds':time.monotonic()-start,'blocks':blocks}
    (r/'comparison.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ('status','cores','mothers','layers','max_death','max_live','elapsed_seconds')}))
if __name__=='__main__':main()
