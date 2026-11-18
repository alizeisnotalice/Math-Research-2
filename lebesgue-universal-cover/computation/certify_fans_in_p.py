#!/usr/bin/env python3
"""Supplement: rigorous fixed-angle proof that the five Sprague fans lie in P.

It imports (silently) the frozen fixed-angle coordinate certificate but does
not modify it.  Each circular arc is covered by rational radial-chord cells.
At an endpoint lying on a P support line, the directed derivative is checked
into the permitted halfplane; all other cells have strict interval slack.
"""
from pathlib import Path
from sigma_parameter import SIGMA_OUT
import contextlib, hashlib, io, json, runpy, sys

here=Path(__file__).parent
with contextlib.redirect_stdout(io.StringIO()):
    d=runpy.run_path(str(here/'certify_sprague_fixed.py'))
I=d['I']; F=d['F']; add=d['add']; sub=d['sub']; mul=d['mul']; dot=d['dot']; norm2=d['norm2']
n=d['n']; m=d['m']
K=d['K']; G=d['G']; Ip=d['Ip']; H=d['H']; J=d['J']; P=d['P']; Qp=d['Qp']; X=d['X']
F3=d['F3']; B3=d['B3']; C3=d['C3']; E3=d['E3']; C2=d['C2']; E2=d['E2']; A1=d['A1']

SUPPORTS=[(f'n{j}',n[j]) for j in range(6)]+[('m3',m[3]),('m5',m[5])]
ARCS=[
    ('CS',F3,K,G,C2,{(0,'m5'):'circle_tangent',(1,'n5'):'transverse'}),
    ('ES-B',B3,Ip,J,E2,{(0,'m3'):'circle_tangent'}),
    ('ES-C',C3,J,H,E2,{(1,'n3'):'circle_tangent'}),
    ('AS-E',E3,Qp,X,A1,{(0,'n1'):'circle_tangent'}),
    ('AS-G',G,X,P,A1,{(1,'n2'):'circle_tangent'}),
]
APEX_ACTIVE={
    'CS': {'n5','m5'}, 'ES-B': {'n3','m3'}, 'ES-C': {'n3','m3'},
    'AS-E': {'n1','n2'}, 'AS-G': {'n1','n2'},
}

def slack(normal,x): return I(F(1,2))-dot(normal,x)
def arc_slack(name,c,u0,u1,tag,normal,active,cells=256):
    """All points c+minor_arc(u0,u1) satisfy normal.x<=1/2."""
    s0=slack(normal,add(c,u0)); s1=slack(normal,add(c,u1)); ss=[s0,s1]
    if active is None:
        assert s0.lo>0 and s1.lo>0, (name,tag,'strict endpoints',s0.out(),s1.out())
    else:
        index,kind=active
        assert ss[index].lo<=0<=ss[index].hi, (name,tag,'specified equality',ss[index].out())
        assert ss[1-index].lo>0, (name,tag,'other endpoint strict',ss[1-index].out())
        if kind=='circle_tangent':
            # The stated projection gives u_index=normal exactly.  Thus
            # normal.(c+u)<=normal.c+1=1/2 over the whole unit circle.
            eq=dot(normal,[u0,u1][index]); assert eq.lo<=1<=eq.hi, (name,tag,'unit tangent identity',eq.out())
            return {'endpoint_slacks':[s0.out(),s1.out()],
                    'mode':'entire-circle tangent inequality',
                    'normal_dot_radius':eq.out()}
    dw=sub(u1,u0); lows=[]; derivatives=[]
    for j in range(cells):
        t=I(F(j,cells),F(j+1,cells)); w=add(u0,mul(t,dw)); r2=norm2(w)
        # du/dt has the following numerator and a positive r^3 denominator.
        du_num=sub(mul(r2,dw),mul(dot(w,dw),w))
        ds_num=-dot(normal,du_num)
        if active==(0,'transverse'):
            # Slack increases monotonically from its specified zero endpoint.
            assert ds_num.lo>0, (name,tag,'inward derivative',j,ds_num.out())
            derivatives.append(ds_num.out()); continue
        if active==(1,'transverse'):
            # Slack decreases monotonically to its specified zero endpoint.
            assert ds_num.hi<0, (name,tag,'inward derivative',j,ds_num.out())
            derivatives.append(ds_num.out()); continue
        u=mul(I(1)/r2.sqrt(),w); s=slack(normal,add(c,u))
        assert s.lo>0, (name,tag,j,s.out()); lows.append(s.lo)
    return {'endpoint_slacks':[s0.out(),s1.out()],
            'cell_slack_lower':str(min(lows)) if lows else None,
            'tangent_derivative_numerators':derivatives}

rows=[]
for name,c,a,b,apex,active_set in ARCS:
    u0,u1=sub(a,c),sub(b,c)
    apex_rows=[]; arc_rows=[]
    for tag,normal in SUPPORTS:
        sa=slack(normal,apex)
        if tag in APEX_ACTIVE[name]:
            assert sa.lo<=0<=sa.hi, (name,tag,'apex equality',sa.out())
        else:
            assert sa.lo>0, (name,tag,'apex strict',sa.out())
        apex_rows.append((tag,sa.out()))
        active=(0,active_set[(0,tag)]) if (0,tag) in active_set else (1,active_set[(1,tag)]) if (1,tag) in active_set else None
        arc_rows.append((tag,arc_slack(name,c,u0,u1,tag,normal,active)))
    rows.append({'fan':name,'apex_slacks':apex_rows,'arc_supports':arc_rows})

out={'status':'passed','scope':'every point of each of the five fixed-angle Sprague fans is in P; convexity of P then covers all fan segments from the apex to its certified arc.  This supplement does not alter or replace the frozen main certificate.',
     'python':sys.version,"sigma_degrees_interval": SIGMA_OUT,'cells_per_arc':256,
     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'base_certificate_sha256':hashlib.sha256((here/'certify_sprague_fixed.py').read_bytes()).hexdigest(),
     'fans':rows}
print(json.dumps(out,indent=2))
