"""Recompute the three preregistered Tomar leaves with the original inputs."""
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
GOAL='.83750'
P_LEAF='180/24225'
FIVE_LEAVES=('1205/140208','1205/140080')
FIVE_NODES=(109641,122279)


def read_rows(path):
    return [json.loads(line) for line in path.read_bytes().splitlines()]


def checked_five(path):
    rows=read_rows(path)
    assert len(rows)==2
    ans={}
    for leaf,nodes,row in zip(FIVE_LEAVES,FIVE_NODES,rows,strict=True):
        assert row['label']==row['leaf']==leaf and row['goal']=='67/80' and row['mode']==4
        assert row['threshold_integer']==120696470013529293
        assert row['stats'][:8]==[1,nodes,(nodes+1)//2,0,0,0,row['stats'][6],0]
        assert row.get('pending',0)==0
        ans[leaf]=dict(complete=True,inputs=5,nodes=nodes,five_terminals=(nodes+1)//2)
    return ans


def checked_six(path,rows):
    result=json.loads(path.read_bytes())
    assert result['leaf']==P_LEAF and result['complete'] and result['pending']==0
    assert result['stitched_runs']==16 and len(rows)==16
    assert result['unresolved_volume_fraction']=='0'
    assert result['nodes']==5615271 and result['terminals']=={
        'five':522841,'P30_cap':2284795,'Q11_cap':0,'mixed':0,
    }
    terminals=sum(result['terminals'].values())
    assert result['nodes']==sum(r['stats'][1] for r in rows)==2*terminals-1
    assert rows[0]['stats'][0:2]==[0,1000000] and rows[0]['pending']==15
    assert all(r['edge_caps'] and r['mode']==0 and r['anchors']==0 and r['guards']==0 and
               r['q_edge_stops']==0 and r['stats'][4]==0 for r in rows)
    assert all(r['stats'][0]==1 and r['pending']==0 for r in rows[1:])
    long_edge=sum(r['p_edge_stops'] for r in rows)
    assert long_edge==1813759
    return dict(complete=True,inputs=6,nodes=result['nodes'],
                five_terminals=result['terminals']['five'],
                P30_terminals=result['terminals']['P30_cap'],P30_long_edge=long_edge,
                exact_box_reassembly=True)


def checked_seven(path,original=None):
    row,=read_rows(path)
    assert row['leaf']==row['label']==P_LEAF and row['mode']==1 and row['goal']=='67/80'
    assert row['threshold_integer']==120696470013529293
    assert row['stats'][:8]==[1,5282119,375679,2212876,52505,0,58,0]
    assert row['stats'][1]==2*sum(row['stats'][2:6])-1
    assert row['pending']==row['anchors']==row['guards']==0 and row['edge_caps']
    assert row['p_edge_stops']==1772333 and row['q_edge_stops']==16008
    if original is not None:
        assert row['lo']==original['lo'] and row['hi']==original['hi']
    return dict(complete=True,leaf=P_LEAF,inputs=7,nodes=row['stats'][1],
                five_terminals=row['stats'][2],P30_terminals=row['stats'][3],
                Q11_terminals=row['stats'][4],P30_long_edge=row['p_edge_stops'],
                Q11_long_edge=row['q_edge_stops'])


def main(a):
    work=a.work.resolve()
    work.mkdir(parents=True,exist_ok=True)
    comp,tree,table,cache=(p.resolve() for p in (a.comp,a.tree,a.table,a.cache))
    python=[sys.executable,'-I','-B']
    digest=lambda p:sha256(p.read_bytes()).hexdigest()
    assert digest(tree/'SHA256SUMS')=='56cbe5dfdadc6b5aa07cf6ec9acf4f637a3c4fa78db54df0ae2d78916d292410'
    # The proof driver verifies the canonical integer arrays, so serialization
    # differences between NumPy versions need not change the mathematical input.
    container_sha=digest(table)

    def run(name,cases,budget,extra):
        frontier=work/(name+'-frontiers')
        argv=[*python,str(HERE/'replay_angle_mesh.py'),str(comp),str(tree),str(table),
              str(work/(name+'-build')),str(cache),'--target',GOAL,'--cases',*cases,'--budget',
              str(budget),'--dense','--frontier-out',str(frontier),*extra]
        output=work/(name+'.jsonl')
        with output.open('w') as raw:
            subprocess.run(argv,stdout=raw,check=True)
        return output,frontier

    def audit(root_record,children,frontiers):
        audit_file=work/'six-audit.json'
        subprocess.run([*python,str(HERE/'audit_partition.py'),str(comp),str(tree),str(audit_file),
                        '--six-only','--receipts',str(root_record),*map(str,children),
                        '--frontier-dirs',*map(str,frontiers)],check=True)
        all_rows=read_rows(root_record)+[x for f in children for x in read_rows(f)]
        return checked_six(audit_file,all_rows)

    folded={}
    seven=None
    if a.audit_frozen:
        src=a.audit_frozen.resolve()
        folded.update(checked_five(src/'five.jsonl'))
        resumes=sorted(src.glob('six-part-?.jsonl'))
        frontier=[src/'six-root-frontiers',*sorted(src.glob('six-part-?-frontiers'))]
        assert len(resumes)==4 and frontier[0].is_dir()
        folded[P_LEAF]=audit(src/'six-root.jsonl',resumes,frontier)
        comparison=src/'seven-control.jsonl'
        if comparison.exists():
            seven=checked_seven(comparison,read_rows(src/'six-root.jsonl')[0])
    else:
        with ThreadPoolExecutor(a.workers) as pool:
            five=None
            if a.stage in {'all','five'}:
                five=pool.submit(run,'five',FIVE_LEAVES,200000,['--modes','4','--panels','0'])
            if a.stage in {'all','six'}:
                root,frontier=run('six-root',[P_LEAF],1000000,
                    ['--panels','32','--q-mesh','--edge-caps','--modes','0','--depth','180',
                     '--reanchor-level','116','--reanchor-stride','12'])
                pending=json.loads((frontier/'180-24225-mode0-frontier.json').read_bytes())
                assert len(pending)==15
                tasks=[]
                cuts=(pending[:9],pending[9:13],pending[13:14],pending[14:])
                for i,chunk in enumerate(cuts):
                    path=work/f'six-part-{i}-input.json'
                    path.write_text(json.dumps(chunk,separators=(',',':'))+'\n')
                    extra=['--subboxes',str(path),'--panels','32','--edge-caps','--modes','0',
                           '--depth','160','--reanchor-level','90','--reanchor-stride','18',
                           '--sweep-limit','800']
                    tasks.append(pool.submit(run,f'six-part-{i}',[P_LEAF],8000000,extra))
                parts,dirs=zip(*(future.result() for future in tasks))
                folded[P_LEAF]=audit(root,parts,[frontier,*dirs])
            if five:
                folded.update(checked_five(five.result()[0]))
            if a.stage=='seven':
                result,_=run('seven-control',[P_LEAF],8000000,
                    ['--panels','32','--q-mesh','--edge-caps','--modes','1','--depth','160',
                     '--reanchor-level','90','--reanchor-stride','18','--sweep-limit','800'])
                seven=checked_seven(result)
    result=dict(scope='three original leaves',target='67/80',areas='strictly greater',
                tree_manifest_sha256=digest(tree/'SHA256SUMS'),table64_container_sha256=container_sha,
                leaves=folded)
    if seven:
        result['seven_control']=seven
    (work/'three-leaf-result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    for name in ('comp','tree','table','work','cache'):
        p.add_argument(name,type=Path)
    p.add_argument('--audit-frozen',type=Path)
    p.add_argument('--stage',choices=('all','five','six','seven'),default='all')
    p.add_argument('--workers',type=int,default=4)
    a=p.parse_args()
    assert a.workers>0
    main(a)
