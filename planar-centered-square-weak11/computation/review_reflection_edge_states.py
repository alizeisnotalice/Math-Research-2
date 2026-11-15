#!/usr/bin/env python3
"""Small exact audit of R46 state meanings and coefficient factors.

This calls the author's kernel only on one-atom period10 toy data, never on
R35/R41. The proof review is independent; full-bbox loops remain tiny.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import hashlib,importlib.util,json,os,platform,sys


def root2_bracket(bits=80):
    D=1<<bits;t=isqrt(2*D*D)
    return F(t,D),F(t+1,D)


def main():
    assert os.path.ismount('/Volumes/data')
    source=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-astra-reflection/planar-centered-square-weak11/computation/reflection_edge_states.py')
    sys.path.insert(0,str(source.parent))
    spec=importlib.util.spec_from_file_location('r46_edge_kernel_audit',source)
    kernel=importlib.util.module_from_spec(spec);spec.loader.exec_module(kernel)
    kernel.INPUT_HASH='one-atom-period10-audit'
    data={'Px':'10','Py':'10','lam':'1','pts':[['0','0','1']]}
    results={}
    exact={'bulk':F(1),'left':F(5,2),'right':F(3,2),'bottom':F(5,2),'top':F(3,2),'seam-old':F(2)}
    seam={'input_sha256':kernel.INPUT_HASH,'rational_reflection_axis':'1/2','c_edge_upper':'1/2'}
    for state in ['bulk','left','right','bottom','top','seam-old','seam-new']:
        result=kernel.calculate(data,state,seam if state.startswith('seam') else None,bits=40)
        lo,hi=F(result['inner_area']),F(result['outer_area'])
        assert lo<=hi
        if state in exact:assert lo<=exact[state]<=hi
        results[state]=result
    rlo,rhi=root2_bracket();glo=F(3,2)-rhi;ghi=F(3,2)-rlo
    def interval(state):return F(results[state]['inner_area']),F(results[state]['outer_area'])
    nlo,nhi=interval('seam-new');olo,ohi=interval('seam-old')
    gamma=(nlo-ohi,nhi-olo)
    assert gamma[0]<=glo<=ghi<=gamma[1]
    Alo,Ahi=interval('bulk')
    Elo=sum(interval(s)[0] for s in ['left','right','bottom','top'])
    Ehi=sum(interval(s)[1] for s in ['left','right','bottom','top'])
    b=(Elo-8*Ahi,Ehi-8*Alo)
    assert b[0]<=0<=b[1]
    assert b[0]+gamma[0]>F(1,20)
    # Move the outward axis by delta. The exact toy answer changes by
    # gamma(c+delta)-gamma(c)=(2-sqrt(2))*delta, not zero.
    delta=F(1,100);moved={'input_sha256':kernel.INPUT_HASH,
       'rational_reflection_axis':str(F(1,2)+delta),'c_edge_upper':'1/2'}
    moved_results={s:kernel.calculate(data,s,moved,bits=40) for s in ['seam-old','seam-new']}
    mnlo,mnhi=map(F,[moved_results['seam-new']['inner_area'],moved_results['seam-new']['outer_area']])
    molo,mohi=map(F,[moved_results['seam-old']['inner_area'],moved_results['seam-old']['outer_area']])
    moved_gamma=(mnlo-mohi,mnhi-molo)
    target_lo=glo+(2-rhi)*delta;target_hi=ghi+(2-rlo)*delta
    assert moved_gamma[0]<=target_lo<=target_hi<=moved_gamma[1]
    assert moved_gamma[0]>gamma[1]
    corner_source=source.with_name('reflection_edge_corners.py')
    cspec=importlib.util.spec_from_file_location('r46_corner_kernel_audit',corner_source)
    corner=importlib.util.module_from_spec(cspec);cspec.loader.exec_module(corner)
    corner.INPUT_HASH=kernel.INPUT_HASH
    corners={}
    expected_corners={'corner-lb':F(25,4),'corner-lt':F(15,4),
                      'corner-rb':F(15,4),'corner-rt':F(9,4),
                      'seam-bottom-old':F(7),'seam-top-old':F(5)}
    for state in ['corner-lb','corner-lt','corner-rb','corner-rt',
                  'seam-bottom-old','seam-bottom-new','seam-top-old','seam-top-new']:
        rr=corner.calculate(data,state,seam if state.startswith('seam') else None,bits=40)
        lo,hi=F(rr['inner_area']),F(rr['outer_area'])
        if state in expected_corners:assert lo<=expected_corners[state]<=hi
        corners[state]=rr
    def ci(state):return F(corners[state]['inner_area']),F(corners[state]['outer_area'])
    C_lo=sum(ci(s)[0] for s in ['corner-lb','corner-lt','corner-rb','corner-rt'])
    C_hi=sum(ci(s)[1] for s in ['corner-lb','corner-lt','corner-rb','corner-rt'])
    assert C_lo<=16<=C_hi
    Gamma_lo=sum(ci('seam-'+s+'-new')[0]-ci('seam-'+s+'-old')[1] for s in ['bottom','top'])
    Gamma_hi=sum(ci('seam-'+s+'-new')[1]-ci('seam-'+s+'-old')[0] for s in ['bottom','top'])
    assert Gamma_lo<=6*glo<=6*ghi<=Gamma_hi
    F6=(2*Elo-32*Ahi+C_lo+Gamma_lo,2*Ehi-32*Alo+C_hi+Gamma_hi)
    assert F6[0]<=6*glo<=6*ghi<=F6[1]
    c0=(16*Alo-4*Ehi+C_lo,16*Ahi-4*Elo+C_hi)
    assert c0[0]<=0<=c0[1]
    out={'status':'PASS','scope':'One-atom toy only; no R35/R41 large state evaluation.',
        'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
        'single_process':True,'bits':40,'reference_sqrt_bits':80,
        'arithmetic':'Fraction and directed dyadic square roots','random_seed':None,'floating_tolerance':None,
        'kernel':str(source),'kernel_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'union_module_sha256':hashlib.sha256(source.with_name('periodic_halo.py').read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'expected_A':'1','expected_E':'8','expected_b':'0','expected_gamma':'3/2-sqrt(2)',
        'gamma_interval':list(map(str,gamma)),'b_interval':list(map(str,b)),
        'outward_axis_delta':str(delta),'expected_moved_gamma':'3/2-sqrt(2)+(2-sqrt(2))/100',
        'moved_gamma_interval':list(map(str,moved_gamma)),
        'states':results,'moved_states':moved_results,'corner_states':corners,
        'corner_kernel_sha256':hashlib.sha256(corner_source.read_bytes()).hexdigest(),
        'expected_Ccorner':'16','expected_Gamma':'6*(3/2-sqrt(2))',
        'Gamma_interval':[str(Gamma_lo),str(Gamma_hi)],'F6_interval':list(map(str,F6))}
    path=Path(__file__).parent/'results/reflection-edge-state-review/exact-toy.json'
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','kernel_sha256':out['kernel_sha256'],'state_calls':17,
                      'max_source_atoms':max(x['window_atoms'] for x in list(results.values())+list(moved_results.values())+list(corners.values())),
                      'max_boxes':max(x['counts']['boxes'] for x in list(results.values())+list(moved_results.values())+list(corners.values())),
                      'output':str(path)},ensure_ascii=False))

if __name__=='__main__':main()
