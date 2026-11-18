"""Replay complete original boxes with piecewise P30 and optional Q11 angle anchors."""
from argparse import ArgumentParser
import ctypes
from fractions import Fraction as F
import json
from math import isqrt
from pathlib import Path
import subprocess
import sys
from time import process_time
import traceback

import numpy as np
from flint import arb, ctx


def support_mesh(proof, lo, hi, panels, fine, area_scale, direction_scale, ceil):
    knots, supports, pads = [], [], []
    for n, axis in ((5, 2), (7, 5), (9, 8)):
        a, b = int(lo[axis]), int(hi[axis])
        assert b - a >= panels
        grid = [a + (b - a) * j // panels for j in range(panels + 1)]
        max_width = max(y - x for x, y in zip(grid[:-1], grid[1:], strict=True))
        pad = ceil(area_scale * F(3, 5) * F(max_width, fine)**2 *
                   F(direction_scale + 2, direction_scale) / 8)
        pads.append(pad)
        knots.extend(grid)
        supports.extend(proof.one(n, F(t, fine)) + pad for t in grid)
    return np.array(knots, np.int64), np.array(supports, np.int64).ravel(), pads


def generate_source(work, base):
    new = base.decode()
    start = new.index("                    L q=p>=h ? hi[k] : lo[k];")
    stop = new.index("\n                }", start)
    old = new[start:stop]
    new = new[:start] + (
        "                    if (mesh.panels) {\n"
        "                        h = mesh.upper(i,j,n,lo[k],hi[k]);\n"
        "                    } else {\n" +
        "\n".join("    "+line for line in old.splitlines()) +
        "\n                    }"
    ) + new[stop:]
    new = new.replace('struct DynamicGuard {',
                      '#include "angle_mesh.hpp"\n'
                      '#ifdef WHOLE_BOX_LONG_EDGE\n#include "long_edge_support.hpp"\n#endif\n'
                      '\nstruct DynamicGuard {', 1)
    new = new.replace("    static constexpr L U=", "    AngleMesh mesh;\n\n    static constexpr L U=", 1)
    new = new.replace('    L collective_boxes=0;', '    L collective_boxes=0,edge_boxes=0;',1)
    start_edge = "        bool needs_collective=false;"
    new = new.replace(start_edge,start_edge+"""
#ifdef WHOLE_BOX_LONG_EDGE
        LongEdgeSupport edge{t};
        std::vector<L> edge_depth;
        bool needs_edge=false;
#endif
""",1)
    end_edge = "            if(!ok)return false;"
    new = new.replace(end_edge,"""
#ifdef WHOLE_BOX_LONG_EDGE
            if(!ok) {
                if(edge_depth.empty()) edge_depth.assign(n,-1);
                for(L k=pointers[i];k<pointers[i+1];++k) {
                    const L* q=triples+8*k;
                    const W raw_gap=W(q[6])*q[7]-W(q[3])*support_upper(q[0])
                        -W(q[4])*support_upper(q[1])-W(q[5])*support_upper(q[2]);
                    if(raw_gap<=0) continue;
                    W raw_sum=0,norm_squares=0;
                    L stops[3]{};
                    for(int j=0;j<3;++j) {
                        const L ind=q[j], high=support_upper(ind);
                        if(edge_depth[ind]<0)
                            edge_depth[ind]=edge.bound(directions[2*ind],directions[2*ind+1],
                                high,e,threshold(ind)-high);
                        const W term=W(q[3+j])*edge_depth[ind];
                        raw_sum+=term;
                        stops[j]=up_signed(term,q[6]);
                        norm_squares+=W(stops[j])*stops[j];
                    }
                    const auto pair=overlap[k];
                    if(pair[0]>=0)norm_squares+=W(2)*stops[pair[0]]*stops[pair[1]];
                    const L gap=L(raw_gap/q[6]);
                    if(raw_gap>raw_sum || (mode>=2 && W(gap)*gap>norm_squares)) {
                        ok=true;
                        needs_edge=true;
                        break;
                    }
                }
            }
#endif
"""+end_edge,1)
    new = new.replace('        if(needs_collective)++collective_boxes;',
                      '        if(needs_collective)++collective_boxes;\n'
                      '#ifdef WHOLE_BOX_LONG_EDGE\n        if(needs_edge)++edge_boxes;\n#endif',1)
    (work/'angle_mesh_dynamic.hpp').write_text(new)


def check_steps(bounds, path):
    box = [list(pair) for pair in bounds]
    assert isinstance(path, list)
    for code in path:
        assert type(code) is int and 0<=code<22
        i, side = divmod(code, 2)
        lo, hi = box[i]
        assert lo<hi and (lo+hi)%2==0
        box[i][1-side] = (lo+hi)//2
    return box


def intersections(guards,tree,index):
    levels=read(index)
    hits=[]
    for part,positions in levels.items():
        for pos in positions:
            ends=[(int(a*2**45),int(b*2**45)) for a,b in tree.leaf(int(part),pos)]
            for row in guards:
                d=[max(a-z,z-b,0) for (a,b),z in zip(ends,row[:-1],strict=True)]
                fits=True
                for j,l in ((0,0),(3,163),(6,115),(9,89)):
                    room=1000*row[-1]-(l*d[j-1] if j else 0)
                    if room<0 or 10**6*(d[j]**2+d[j+1]**2)>room*room:
                        fits=False
                        break
                if fits:
                    hits.append([int(part),pos])
                    break
    return dict(overlap=len(hits),all_leaves=sum(map(len,levels.values())),leaves=hits)


def main(a):
    ctx.prec=160
    sys.path.insert(0,str(a.comp))
    from build_q11_leaf_bank import reflected_average, EPS
    from replay_q11_with_p30 import PointProof
    from certify_tomar_six_high import FINE, Q, U, SourceTree, build_kernel, ceil, digest, load_table
    from certify_nine import ar

    source=(a.comp/'tomar_six_083730_guard.cpp').read_bytes()
    assert digest(source)=="42a4ee583d22940b3cd0b05f04285ba7e3d351745de707476211dce4efba5f81"
    pi=arb.pi()
    r=ar(F(".707"))
    assert (pi-(1/(2*r)).acos())/4+(r*r-ar(F(1,4))).sqrt()/2 > ar(a.target)
    assert r*r < ar(F(1,2))
    data=read(a.comp/'tomar-six-high-input.json')
    tree=SourceTree(a.tree)
    table=load_table(a.table,tree,data['tables']['64'])
    pcloud=read(a.comp/'finite-point-leaf-input.json')
    pcloud['cells']=read(a.comp/'finite-point-global-input.json')['cells']
    pp=PointProof(pcloud)
    if any(mode&1 for mode in a.modes):
        qcloud,qden=reflected_average(read(a.comp/'multi-host-points.json'))
        bank=read(a.comp/'q11-leaf-bank.json')
        qq=PointProof({**bank,'points':qcloud,'denominator':qden,'direction_scale':U},EPS)
    else:
        qq=pp
    for n,c in ((5,163),(7,115),(9,89)):
        assert (pi/(2*n)).tan()/2 < ar(F(c,1000))
    work=a.work.resolve()
    build_kernel(work,data,a.cache)
    lib=work/'guard.so'
    generate_source(work,source)
    extra=Path(__file__).resolve().with_name('angle_mesh_leaf_guard.cpp')
    edge_flags=['-DWHOLE_BOX_LONG_EDGE'] if a.edge_caps else []
    subprocess.run(['clang++','-O3','-std=c++17','-ftrapv','-fPIC','-shared',*edge_flags,
                    f'-fmodules-cache-path={a.cache}','-I',str(work),'-I',str(a.comp),
                    '-I',str(Path(__file__).resolve().parent),
                    str(extra),'-o',str(lib)],check=True)
    call=ctypes.CDLL(str(lib)).mixed_point_leaf_run
    vec=lambda d:np.ctypeslib.ndpointer(dtype=d,ndim=1,flags='C_CONTIGUOUS')
    i64,i32,M=vec(np.int64),vec(np.int32),ctypes.c_int
    pointer=ctypes.POINTER(ctypes.c_int64)
    anchor_callback=ctypes.CFUNCTYPE(M,pointer,pointer,pointer,pointer)
    center_callback=ctypes.CFUNCTYPE(ctypes.c_int64,pointer,pointer)
    call.argtypes=[i64,i64,i64,i32,i32,ctypes.c_int64,M,M,M,*[i64]*5,M,*[i64]*5,M,i64,M,i64,M,M,M,i64,
                   i64,i64,i64,M,M,i64,anchor_callback,M,M,center_callback,M,M,M]
    call.restype=M

    target=ceil(a.target*2*Q*Q)
    assert F(target,2*Q*Q)>a.target
    def scaled(vals):
        rationals=[F(z)*FINE for z in vals]
        assert all(z.denominator==1 for z in rationals)
        return np.array(list(map(int,rationals)),dtype=np.int64)
    original=scaled(tree.root['hi'])
    before=process_time()
    if a.cert is not None:
        assert all(mode&1 for mode in a.modes), "Q11 center proofs require the seven-input mode"
        proposed=read(a.cert)
        if proposed.get('angle_verifier')=='events':
            sys.path.insert(0,str(Path(__file__).resolve().parent))
            from compress_q11_event_banks import replay
            certs=replay(a.comp,proposed)
        else:
            from certify_mixed_leaf_centers import replay_bank
            certs=replay_bank(a.comp,proposed)
    else:
        certs=None
    bank_cpu=process_time()-before
    if certs:
        assert F(certs['target'])==a.target and certs['scale']==2*FINE
        assert certs['Q11_sha256']==digest((a.comp/'multi-host-points.json').read_bytes())
    p_certs=None
    if a.p_bank:
        sys.path.insert(0,str(Path(__file__).resolve().parent))
        from p30_center import replay as replay_p_bank
        p_certs=replay_p_bank(a.comp,read(a.p_bank))
        assert F(p_certs['target'])==a.target
    cut_banks={}
    if a.cut_bank:
        sys.path.insert(0,str(Path(__file__).resolve().parent))
        from p30_translation_cuts import replay as replay_cut
        for cut_source in a.cut_bank:
            data_cut=read(cut_source)
            assert F(data_cut['target'])==a.target
            center_row,_=replay_cut(a.comp,data_cut)
            cut_banks.setdefault(data_cut['leaf'],[]).append(center_row+[0])
    if a.center_limit:
        sys.path.insert(0,str(Path(__file__).resolve().parent))
        from p30_center import load_tools, make_one, replay as replay_p_bank
        center_tools=load_tools(a.comp)
        assert a.target==F('.83750')
        a.dynamic_center_out.mkdir(parents=True,exist_ok=True)
    if a.overlap_index is not None:
        assert certs and a.reuse
        print(json.dumps(intersections(certs['leaves'][a.reuse]['guards'],tree,a.overlap_index)),flush=True)
    requested=[(key,key,None,[]) for key in a.cases]
    if a.subboxes:
        requested=[]
        for item in read(a.subboxes):
            assert item['leaf'] in a.cases
            requested.append((item['label'],item['leaf'],item['bounds'],item.get('steps',[])))
    if a.frontier_out:
        a.frontier_out.mkdir(parents=True,exist_ok=True)
    for label,key,restriction,previous_step in requested:
        box=tree.leaf(*map(int,key.split('/')))
        root_bounds=[[int(l*FINE),int(h*FINE)] for l,h in box]
        if restriction:
            clipped=[tuple(F(v,FINE) for v in pair) for pair in restriction]
            assert all(l<=a0<b0<=h for (l,h),(a0,b0) in zip(box,clipped,strict=True))
            if previous_step:
                assert check_steps(root_bounds,previous_step)==restriction
            box=clipped
        lo,hi=(scaled(t) for t in zip(*box,strict=True))
        before=process_time()
        d1,d2=pp.presolve(box),qq.presolve(box)
        current_table=table
        if a.dense:
            from tomar_p30_scaled_material import scaled_material
            current_table=scaled_material(table,lo,hi,original)
        knots=p_values=q_values=np.empty(0,np.int64)
        pads=[]
        if a.panels:
            knots,p_values,pads=support_mesh(pp,lo,hi,a.panels,FINE,Q,U,ceil)
            if a.q_mesh and any(mode&1 for mode in a.modes):
                q_knots,q_values,q_pads=support_mesh(qq,lo,hi,a.panels,FINE,Q,U,ceil)
                assert np.array_equal(knots,q_knots) and pads==q_pads
        prep_cpu=process_time()-before
        records=certs['leaves'].get(a.reuse or key,{}) if certs else {}
        material=[r+[int(c.get('strategy','mixed')=='cap')]
                  for r,c in zip(records.get('guards',[]),records.get('certs',[]),strict=True)]
        if p_certs and key in p_certs['leaves']:
            material.extend(r+[0] for r in p_certs['leaves'][key]['guards'])
        material.extend(cut_banks.get(key,[]))
        coords=np.array(material,np.int64).ravel()
        assert len(coords)%13==0
        for mode in a.modes:
            if mode&2:
                assert len(coords)>0
            out=np.zeros(42+len(material),np.int64)
            trace=np.zeros(22*a.trace,np.int64)
            pending=np.zeros((a.depth+max(a.sweep_limit,1)+1,a.depth+23),np.int64)
            new_lot=dict(guards=[],certs=[])
            new_bank_path=None
            if a.center_limit:
                new_bank=dict(target=str(a.target),scale=2*FINE,
                              P30_sha256=center_tools[-1],leaves={key:new_lot})
                new_bank_path=a.dynamic_center_out/(label.replace('/','-')+f'-mode{mode}-p30-bank.json')

            @anchor_callback
            def reanchor(a_lo,a_hi,p_out,q_out):
                try:
                    new_box=[(F(a_lo[k],FINE),F(a_hi[k],FINE)) for k in range(11)]
                    p_ends=pp.presolve(new_box)[2]
                    np.ctypeslib.as_array(p_out,shape=p_ends.shape)[:]=p_ends
                    if mode&1:
                        q_ends=qq.presolve(new_box)[2]
                        np.ctypeslib.as_array(q_out,shape=q_ends.shape)[:]=q_ends
                    return 0
                except Exception:
                    traceback.print_exc()
                    return 1

            @center_callback
            def make_center(a_lo,a_hi):
                try:
                    center=[int(a_lo[k]+a_hi[k]) for k in range(11)]
                    d=[int(a_hi[k]-a_lo[k]) for k in range(11)]
                    physical=[]
                    for k,slope in ((0,0),(3,163),(6,115),(9,89)):
                        trans=isqrt(d[k]**2+d[k+1]**2)+1
                        physical.append(trans+(slope*d[k-1]+999)//1000)
                    radius=F(max(physical)+1,2*FINE)
                    high=max(F('.00012'),radius*2)
                    choices=[high,max(F('.00006'),radius*F(3,2)),radius]
                    name=f'{label}-generated-{len(new_lot["guards"])}'
                    answer=make_one(center_tools,center,list(dict.fromkeys(choices)),name)
                    if answer is None:
                        return 0
                    checked_row,proof=answer
                    new_lot['guards'].append(checked_row)
                    new_lot['certs'].append(proof)
                    new_bank_path.write_text(json.dumps(new_bank,separators=(',',':'))+'\n')
                    print(f'{name} radius={float(F(checked_row[-1], 2*FINE)):.9g}',
                          file=sys.stderr,flush=True)
                    return checked_row[-1]
                except Exception:
                    traceback.print_exc()
                    return -1

            before=process_time()
            assert call(lo,hi,original,*current_table,target,a.depth,a.budget,mode,
                        *d1,len(d1[1]),*d2,len(d2[1]),coords,len(coords)//13,
                        trace,a.trace,a.stride,a.trace_kind,out,
                        knots,p_values,q_values,a.panels,int(a.q_mesh),
                        pending.ravel(),reanchor,a.reanchor_level,a.reanchor_stride,
                        make_center,a.center_level,a.center_limit,a.sweep_limit)==0
            used=process_time()-before
            stats=out[:35].tolist()
            uses=out[35:35+len(material)].tolist()
            pending_count,anchors,center_calls,center_stops,p_edge,q_edge,pending_visited=map(int,out[-7:])
            assert stats[5]==sum(stats[33:35])==sum(uses)+center_stops
            assert len(new_lot['guards'])<=center_calls<=a.center_limit
            if new_bank_path and new_lot['guards']:
                replay_p_bank(a.comp,read(new_bank_path))
            terminal_count=sum(stats[2:6])
            if stats[0]:
                assert stats[1]==2*terminal_count-1 and stats[7]==0 and pending_count==0
            else:
                assert 1<=pending_count<=len(pending) and 1<=stats[7]<=6
                assert 0<=pending_visited<=pending_count
                assert stats[1]==2*terminal_count+pending_count-1+pending_visited
            suffixes=[]
            frontier=[]
            bounds=np.column_stack((lo,hi)).tolist()
            for i,z in enumerate(pending[:pending_count]):
                length=int(z[0])
                suffix=z[23:23+length].tolist()
                end=np.column_stack((z[1:12],z[12:23])).tolist()
                assert check_steps(bounds,suffix)==end
                suffixes.append(tuple(suffix))
                address=previous_step+suffix
                frontier.append(dict(leaf=key,label=f'{label}-{i}',bounds=end,steps=address))
            assert len(set(suffixes))==len(suffixes)
            sorted_suffixes=sorted(suffixes)
            assert all(x[:len(y)]!=y for y,x in zip(sorted_suffixes,sorted_suffixes[1:]))
            if pending_count:
                assert np.array_equal(pending[0,1:23],stats[11:33])
            if not mode&1:
                assert stats[4]==0
            if not mode&2:
                assert stats[5]==0
            row=dict(leaf=key,label=label,goal=str(a.target),threshold_integer=target,
                     lo=lo.tolist(),hi=hi.tolist(),panels=a.panels,pads=pads,q_mesh=a.q_mesh,
                     pending=pending_count,anchors=anchors,
                     pending_visited=pending_visited,sweep_limit=a.sweep_limit,
                     reanchor_level=a.reanchor_level,reanchor_stride=a.reanchor_stride,
                     center_level=a.center_level,center_limit=a.center_limit,
                     center_calls=center_calls,center_stops=center_stops,generated=len(new_lot['guards']),
                     cut_guards=len(cut_banks.get(key,[])),
                     edge_caps=a.edge_caps,p_edge_stops=p_edge,q_edge_stops=q_edge,
                     mode=mode,budget=a.budget,stats=stats,dense=a.dense,
                     prep_cpu=prep_cpu,search_cpu=used,bank_cpu=bank_cpu,guards=len(coords)//13,guard_uses=uses)
            print(json.dumps(row),flush=True)
            if a.frontier_out:
                path=a.frontier_out/(label.replace('/','-')+f'-mode{mode}-frontier.json')
                path.write_text(json.dumps(frontier,separators=(',',':'))+'\n')
            if a.trace and (mode&1 or a.trace_kind==1):
                a.work.joinpath(key.replace('/','-')+f'-{mode}-trace.json').write_text(
                    json.dumps(trace[:22*stats[10]].reshape(-1,22).tolist(),separators=(',',':'))+'\n')


def read(path):
    return json.loads(path.read_bytes())


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    for key in ('comp','tree','table','work','cache'):
        p.add_argument(key,type=Path)
    p.add_argument('--cert',type=Path)
    p.add_argument('--p-bank',type=Path)
    p.add_argument('--cut-bank',nargs='+',type=Path)
    p.add_argument('--overlap-index',type=Path)
    p.add_argument('--reuse',help='reuse another leaf\'s certified geometric guards')
    p.add_argument('--target',type=F,required=True)
    p.add_argument('--budget',type=int,required=True)
    p.add_argument('--depth',type=int,default=116)
    p.add_argument('--modes',type=int,choices=(0,1,2,3,4),nargs='+',required=True)
    p.add_argument('--trace',type=int,default=0)
    p.add_argument('--trace-kind',type=int,choices=(0,1),default=0)
    p.add_argument('--stride',type=int,default=1)
    p.add_argument('--subboxes',type=Path)
    p.add_argument('--frontier-out',type=Path)
    p.add_argument('--dense',action='store_true')
    p.add_argument('--panels',type=int,default=32)
    p.add_argument('--q-mesh',action='store_true')
    p.add_argument('--edge-caps',action='store_true')
    p.add_argument('--sweep-limit',type=int,default=0)
    p.add_argument('--reanchor-level',type=int,default=0)
    p.add_argument('--reanchor-stride',type=int,default=16)
    p.add_argument('--center-level',type=int,default=128)
    p.add_argument('--center-limit',type=int,default=0)
    p.add_argument('--dynamic-center-out',type=Path)
    p.add_argument('--cases',nargs='+',required=True)
    args=p.parse_args()
    assert args.stride>0 and args.trace>=0 and args.panels>=0
    assert 0<=args.reanchor_level<=args.depth and args.reanchor_stride>0
    assert args.center_level>=0 and args.center_limit>=0
    assert args.sweep_limit>=0
    assert bool(args.center_limit)==bool(args.dynamic_center_out)
    main(args)
