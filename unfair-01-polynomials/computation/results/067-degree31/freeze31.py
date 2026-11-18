"""Bind complete original and checked degree-31 ledgers and freeze part files."""
from argparse import ArgumentParser
from collections import Counter
import gzip
from hashlib import sha256
import json
from pathlib import Path

from census31 import OLD, sha


def same_receipts(original,checked,candidates):
    with original.open() as a, checked.open() as b:
        start=json.loads(next(a))
        confirm=json.loads(next(b))
        identity=('start','stop','census_sha256','sources_sha256')
        assert all(start[k]==confirm[k] for k in identity)
        assert start['mode']=='generate' and confirm['mode']=='verify'
        assert confirm['witness_sha256']==sha(original)
        counts=Counter()
        unresolved=[]
        for n in range(start['start'],start['stop']):
            previous,independent=next(a),next(b)
            assert previous==independent
            row=json.loads(previous)
            assert row['index']==n and row['shape']==candidates[n]['shape']
            counts[row['status']]+=1
            if row['status']=='unit':
                counts['mod'+str(row['prime'])]+=1
                counts[row['method']]+=1
            else:
                assert row['status']=='unresolved'
                unresolved.append(row)
        end,other=json.loads(next(a)),json.loads(next(b))
        assert end['status']==other['status']=='complete' and end['counts']==other['counts']==dict(counts)
        assert end['start']==other['start']==start['start'] and end['stop']==other['stop']==start['stop']
        assert not a.read(1) and not b.read(1)
    return start,end,other,counts,unresolved


def archive(original,target):
    files=[]
    unpacked=sha256()
    with original.open('rb') as f:
        i=0
        while raw:=f.read(6_000_000):
            name=f'{original.stem}-{i:02d}.part.gz'
            encoded=gzip.compress(raw,compresslevel=6,mtime=0)
            assert len(encoded)<10*1024**2
            with (target/name).open('xb') as part:
                part.write(encoded)
            data=gzip.decompress((target/name).read_bytes())
            assert data==raw
            unpacked.update(data)
            files.append({'file':name,'raw_bytes':len(raw),'compressed_bytes':len(encoded),
                          'raw_sha256':sha256(raw).hexdigest(),'compressed_sha256':sha256(encoded).hexdigest()})
            i+=1
    assert unpacked.hexdigest()==sha(original)
    return files


def check_real(row,real):
    if row['shape']!=real['shape']:
        return False
    assert row['exact_system']['variable_indices']==[8,9]
    terms=dict(row['exact_system']['early_zero_equations'])[23]
    known=dict((tuple(mon),coef) for mon,coef in terms)
    expected=list(zip(((2,0),(1,1),(1,0),(0,2),(0,1),(0,0)),real['row23_coefficients'],strict=True))
    assert known==dict(expected) and real['strict_lower_bound']=='11/160'
    return True


def main():
    if not __debug__:
        raise SystemExit('assertions required')
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--census',type=Path,required=True)
    parser.add_argument('--census-audit',type=Path,required=True)
    parser.add_argument('--real',type=Path,required=True)
    parser.add_argument('--units',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args()
    original=json.loads(a.census.read_bytes())
    audit=json.loads(a.census_audit.read_bytes())
    real=json.loads(a.real.read_bytes())
    assert real['status']==audit['status']=='verified'
    assert audit['census_sha256']==sha(a.census)
    labs=original['ledger']['candidates']
    unitfiles=sorted(a.units.glob('raw-*.jsonl'))
    assert unitfiles
    a.out.mkdir(exist_ok=False)
    accepted=Counter()
    slices=[]
    unresolved=[]
    resolved=[]
    source={}
    end=0
    for raw in unitfiles:
        checked=raw.with_name(raw.name.replace('raw-','check-',1))
        header,terminal,replay,counts,unknowns=same_receipts(raw,checked,labs)
        assert header['start']==end and header['census_sha256']==sha(a.census)
        end=header['stop']
        accepted.update(counts)
        for name,digest in header['sources_sha256'].items():
            assert source.setdefault(name,digest)==digest
        for row in unknowns:
            if check_real(row,real):
                resolved.append({'shape':row['shape'],'index':row['index'],'strictly_positive_residual':real['strict_lower_bound']})
            else:
                unresolved.append(row)
        slices.append({'start':header['start'],'stop':header['stop'],'counts':terminal['counts'],
                       'generation_seconds':terminal['seconds'],'identity_replay_seconds':replay['seconds'],
                       'raw_sha256':sha(raw),'identity_sha256':sha(checked),'parts':archive(raw,a.out)})
    assert end==len(labs)==accepted['unit']+accepted['unresolved']
    assert accepted['unresolved']==len(resolved)+len(unresolved)
    for supplement in (audit,real):
        for name,digest in supplement['sources_sha256'].items():
            assert source.setdefault(name,digest)==digest
    here=Path(__file__).resolve().parent
    for name,digest in source.items():
        location=here/name if (here/name).exists() else OLD/name
        assert sha(location)==digest
    count=original['counts']
    assert count['all']==count['early_rejections']+count['late_rejections']+accepted['unit']+len(resolved)+len(unresolved)
    result={'status':'verified' if not unresolved else 'incomplete','degree':31,
            'counts':{'shapes':count,'modular':dict(accepted),'real_certificate':len(resolved),'unresolved_final':len(unresolved)},
            'census':{'sha256':sha(a.census),'parts':archive(a.census,a.out),'seconds':original['seconds']},
            'census_audit':{'sha256':sha(a.census_audit),'seconds':audit['seconds']},
            'real_check':{'sha256':sha(a.real),'resolved':resolved},'unresolved':unresolved,
            'slices':slices,'python':header['python'],'sympy':header['sympy'],'sources_sha256':source}
    with (a.out/'verification.json').open('x') as f:
        f.write(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],result['counts'],flush=True)


if __name__=='__main__':
    main()
