"""Replay the two R228 area and pose results using the public project geometry."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace

PACKAGE=Path(__file__).resolve().parent
INPUT=PACKAGE/'inputs'

SINGLE={
    'single-positive':('first-positive-19-strong.json',300,23,100000),
    'single-reflected':('reflect-pair-19.json',350,24,170000),
    'single-both':('pair-19.json',300,24,100000),
}
PLACED={
    'single-five':'reflection-escape-five-19.json',
    'wrong-four':'negative-four.json',
    'wrong-plus-Y':'repaired-negative.json',
    'reseated-five':'positive-five.json',
    'free-regular-five':'regular-five.json',
}
STAGES=('base-host','single-positive','single-reflected','single-five',
        'wrong-four','wrong-plus-Y','reseated-five','free-regular-five',
        'free-fixed-lower','free-B-lower','wrong-material-lower','free-angle-box','single-both')
B=F('.844082046521218581802948589770767704354')

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def read(path):
    return json.loads(path.read_bytes())

def write(path, value):
    path.write_text(json.dumps(value,indent=2)+'\n')

def run_stage(comp, work, out, stage):
    rel=comp/'released-hexagon'
    repair=comp/'single-repair'
    host=rel/'single-H1-input.json'
    sys.path[:0]=[str(PACKAGE),str(rel),str(comp)]
    if stage=='base-host':
        from compare import run
        result=run(comp,repair,rel/'continuous-153-input.json',host)
        write(out,result)
    elif stage in SINGLE:
        from certify_joint import run
        source,steps,depth,budget=SINGLE[stage]
        run(SimpleNamespace(comp=comp,host=host,fixture=INPUT/'joint-fixture.json',
            pair=INPUT/source,single_repair=repair,output=out,steps=steps,depth=depth,budget=budget))
    elif stage in PLACED:
        from certify_hull import prove
        write(out,prove(comp,host,INPUT/PLACED[stage],repair,rel))
    elif stage in ('free-fixed-lower','free-B-lower'):
        from mixed_level import run
        fixed=stage=='free-fixed-lower'
        write(out,run(SimpleNamespace(comp=comp,single=rel,difference=work/'base-host.json',
            area_target=F('.844082054') if fixed else B,cap_replay=work/'single-positive.json',
            negative=F('3.1199e-7') if fixed else F('2.8e-7'),steps=180 if fixed else 200,
            free_depth=6 if fixed else 7,witness=comp/'free-face-expansion'/'free-witness.json')))
    elif stage=='wrong-material-lower':
        from material_lower import run
        write(out,run(comp,rel,INPUT/'negative-four.json',work/'free-fixed-lower.json',work/'base-host.json'))
    elif stage=='free-angle-box':
        from phase_box import run
        write(out,run(comp,rel,work/'free-B-lower.json'))

def row_bodies(recipe):
    return [(F(row[0]),F(row[1])) for row in recipe['test_bodies']]

def validate_single(work, lo):
    from flint import arb
    p=read(work/'single-positive.json')
    r=read(work/'single-reflected.json')
    u=read(work/'single-five.json')
    for result,sign,gain in ((p,1,F('1.30e-7')),(r,-1,F('9.22e-8'))):
        assert result['complete'] and result['pending_states']==0 and result['unresolved']==[]
        assert result['first_handedness']==sign and F(result['relative_gain'])==gain
        assert result['state_counts']['jointAH']>0
        assert result['cap_area_upper']==p['cap_area_upper']
    lower=F('.8440820130770037')
    upper=F('.8440820140975347851512')
    assert F('.8440819208770037')<lo
    assert F('.8440819208770037')+F(r['relative_gain'])==lower
    assert lo+F(p['relative_gain'])>B
    assert lower<F(u['area_upper'])<upper<B
    assert u['complete'] and list(map(F,u['facet_budget']))==[0,F('3.01e-7'),0,0,0,0]
    assert all(arb(x)<0 for row in u['oblique_excess'] for x in row)
    recipe=read(INPUT/PLACED['single-five'])
    assert len(recipe['test_bodies'])==5 and all(row[2:]==['239.999999',-1] for row in recipe['test_bodies'])
    ids=row_bodies(recipe)
    assert row_bodies(read(INPUT/SINGLE['single-positive'][0]))==[ids[i] for i in (0,2,3)]
    assert row_bodies(read(INPUT/SINGLE['single-reflected'][0]))==[ids[i] for i in (0,4)]
    assert row_bodies(read(INPUT/SINGLE['single-both'][0]))==ids[:2]
    return dict(scope='D contains C* and keeps the four pinned faces and H4, with H1 released by 3.01e-7',
        lower=str(lower),strict_upper=str(upper),true_reflected_upper=u['area_upper'],
        B_minus_upper=str(B-upper),B_minus_lower=str(B-lower),
        first_positive_gain=p['relative_gain'],first_reflected_gain=r['relative_gain'],
        joint_leaves=[p['state_counts']['jointAH'],r['state_counts']['jointAH']])

def validate_free(work):
    from flint import arb
    f=read(work/'free-fixed-lower.json')
    b=read(work/'free-B-lower.json')
    phase=read(work/'free-angle-box.json')
    material=read(work/'wrong-material-lower.json')
    u={name:read(work/(name+'.json')) for name in PLACED if name!='single-five'}
    assert F(f['area_budget'])==F('.844082054')>B and F(f['forbidden_negative_magnitude'])==F('3.1199e-7')
    assert F(b['area_budget'])==B and F(b['forbidden_negative_magnitude'])==F('2.8e-7')
    assert f['witness']==b['witness']==phase['extra_body']
    assert F(f['cap_upper'])==F(b['cap_upper'])
    assert arb(u['wrong-four']['width_excess'][0])<-arb('3.1199e-7')
    assert F('.8440820083564575')<F(material['area_lower'])<F(u['wrong-four']['area_upper'])<F('.8440820090895462')<B
    assert F('.844082054')<F(u['wrong-plus-Y']['area_upper'])<F('.8440820774558181')
    assert F(u['reseated-five']['area_upper'])<F('.8440820067607142')<B
    assert F(u['free-regular-five']['area_upper'])<F('.8440820197824776')<B
    assert all(z['complete'] and all(arb(gap)<0 for row in z['oblique_excess'] for gap in row) for z in u.values())
    expected={'wrong-four':[0,F('3.01e-7'),0,0,F('3.21e-7'),0],
              'wrong-plus-Y':[0,F('3.01e-7'),0,0,F('3.21e-7'),0],
              'reseated-five':[0,F('3.01e-7'),0,F('2e-8'),0,0],
              'free-regular-five':[0,F('3.01e-7'),0,0,0,0]}
    assert all(list(map(F,u[k]['facet_budget']))==x for k,x in expected.items())
    files={name:read(INPUT/path) for name,path in PLACED.items()}
    wrong=files['wrong-four']['test_bodies']
    assert files['wrong-plus-Y']['test_bodies'][:4]==wrong
    physical=row_bodies(files['wrong-plus-Y'])
    assert len(physical)==5 and row_bodies(files['wrong-four'])==physical[:4]
    assert all(row_bodies(files[name])==physical for name in ('reseated-five','free-regular-five'))
    single_bodies=row_bodies(files['single-five'])
    assert len(single_bodies)==5 and single_bodies[:4]==physical[:4] and single_bodies[-1]!=physical[-1]
    assert physical[-1]==(F(f['witness']['w']),F(f['witness']['e']))
    assert phase['all_angles'] and phase['regions_degrees']=={
        '1':[['0','.00085'],['354.01','354.28'],['359.9975','360']],
        '-1':[['239.9975','240.00085']],
    }
    assert list(map(F,phase['parameter_radius']))==[F('4e-7'),F('1e-7')]
    assert F(phase['q_negative_lower'])==F('2.8e-7')
    for z in (f,b):
        assert set(z['angle_cover'])=={'1','-1'}
        assert all(i['leaves']['joint_translation_partition']>0 for i in z['angle_cover'].values())
    return dict(scope='fixed D_- for free repair lower; only C* retained in contrasting reseating; Y differs from single-face P5',
        wrong_four_lower=material['area_lower'],wrong_four_upper=u['wrong-four']['area_upper'],
        fixed_repair_lower=f['area_budget'],fixed_repair_upper=u['wrong-plus-Y']['area_upper'],
        reseated_upper=u['reseated-five']['area_upper'],
        single_face_reseated_upper=u['free-regular-five']['area_upper'],
        width_bounds={i:z['facet_budget'] for i,z in u.items()},
        q_strict_lower=str(-F(b['forbidden_negative_magnitude'])),
        allowed_angle_intervals_degrees=phase['regions_degrees'],angle_box_halfwidth=phase['parameter_radius'],
        folded_first_q=u['wrong-four']['width_excess'][0])

def orchestrate(comp, work, out, section, jobs):
    from flint import __version__ as flint_version
    assert flint_version=='0.9.0',flint_version
    work.mkdir(parents=True,exist_ok=True)
    run_dir=Path(tempfile.mkdtemp(prefix='r228-',dir=work))
    logs={}
    started=time.monotonic()

    def launch(stage):
        begin=time.monotonic()
        receipt=run_dir/(stage+'.json')
        log=run_dir/(stage+'.log')
        args=[sys.executable,'-I','-B',str(Path(__file__).resolve()),str(comp),str(run_dir),str(receipt),'--stage',stage]
        print('start',stage,flush=True)
        with log.open('w') as dest:
            subprocess.run(args,stdout=dest,stderr=subprocess.STDOUT,check=True)
        print('verified',stage,round(time.monotonic()-begin,2),flush=True)
        return stage,dict(receipt=receipt.name,sha256=digest(receipt),stdout=log.name,
                          seconds=round(time.monotonic()-begin,4))

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        def wave(names):
            for stage,result in pool.map(launch,names):logs[stage]=result
        wave(('base-host',))
        wave(('single-positive','single-reflected','single-five') if section in ('all','single') else ('single-positive',))
        if section in ('all','free'):
            wave(('free-fixed-lower','free-B-lower','wrong-four','wrong-plus-Y','reseated-five','free-regular-five'))
            wave(('wrong-material-lower','free-angle-box'))
        if section=='pair':
            wave(('single-both',))

    baseline=read(run_dir/'base-host.json')
    lo=F(baseline['single']['lower'])
    hi=F(baseline['single']['upper'])
    assert F('.8440819208770037')<lo<hi<F('.8440819212529031')
    if section in ('all','single'):
        single=validate_single(run_dir,lo)
    else:single=None
    free=validate_free(run_dir) if section in ('all','free') else None
    if section=='pair':
        re=read(run_dir/'single-both.json')
        assert re['complete'] and not re['pending_states'] and not re['unresolved'] and re['first_handedness'] is None
        assert F(re['relative_gain'])==F('8.89e-8') and re['state_counts']['jointAH']>0
    paths=[*PACKAGE.glob('*.py'),*INPUT.glob('*.json')]
    frozen={str(i.relative_to(PACKAGE)):digest(i) for i in sorted(paths)}
    write(out,dict(complete=True,mathematical_replay=True,section=section,python=platform.python_version(),
        flint=flint_version,arb_bits=192,run_directory=str(run_dir),elapsed_seconds=round(time.monotonic()-started,3),
        executed=logs,package_sources_and_inputs=frozen,base_host_area_exact=[str(lo),str(hi)],
        single_face=single,free_expansion=free))
    print('completed',section,str(out),flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('comp','work','out'):p.add_argument(name,type=Path)
    p.add_argument('--section',choices=('all','single','free','pair'),default='all')
    p.add_argument('--jobs',choices=range(1,5),type=int,default=2)
    p.add_argument('--stage',choices=STAGES,help=argparse.SUPPRESS)
    a=p.parse_args()
    comp,work,out=(q.resolve() for q in (a.comp,a.work,a.out))
    assert not comp.is_relative_to(PACKAGE)
    assert not work.is_relative_to(comp) and not out.is_relative_to(comp), 'Public computation directory is read-only.'
    assert __debug__, 'Run with Python assertions enabled.'
    if a.stage:
        run_stage(comp,work,out,a.stage)
    else:
        orchestrate(comp,work,out,a.section,a.jobs)

if __name__=='__main__':
    main()
