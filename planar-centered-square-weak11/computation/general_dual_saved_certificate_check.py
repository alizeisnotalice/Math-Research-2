"""Independent exact reread of archived P=M certificates, without LP search.

For each positive coefficient the equality certificate forces c_S=m(S).
To verify this geometrically, find one radius in the continuous arrangement
layer immediately below sqrt(m(S))/2. Exact inclusion-exclusion computes
positive center area that captures S and excludes every other atom.
Constancy of incidence on that layer proves c_S=m(S). This does not call the
probe's coefficient enumeration or any optimization library.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
from datetime import datetime,timezone
import hashlib,itertools,json,math,platform,sys,time


def capacity_saturation_check(points,masses,S):
    mass=sum(masses[i] for i in S)
    critical={F(0)}
    for p,q in itertools.combinations(points,2):
        for axis in (0,1):
            r=abs(p[axis]-q[axis])/2
            if 4*r*r<mass:
                critical.add(r)
    last=max(critical)
    bits=12
    while True:
        D=1<<bits
        q=math.isqrt((mass.numerator*D*D)//(4*mass.denominator))
        radius=F(q-1,D)
        if radius>last:
            break
        bits+=8
        assert bits<1024
    assert 0<radius and 4*radius*radius<mass
    # Feasible centers for including S. Excluding j removes its closed box.
    base=(max(points[i][0] for i in S)-radius,
          min(points[i][0] for i in S)+radius,
          max(points[i][1] for i in S)-radius,
          min(points[i][1] for i in S)+radius)
    def intersect_area(boxes):
        x0=max(b[0] for b in boxes);x1=min(b[1] for b in boxes)
        y0=max(b[2] for b in boxes);y1=min(b[3] for b in boxes)
        return max(F(0),x1-x0)*max(F(0),y1-y0)
    outside=[i for i in range(len(points)) if i not in S]
    excluded={i:(points[i][0]-radius,points[i][0]+radius,
                 points[i][1]-radius,points[i][1]+radius) for i in outside}
    available=intersect_area([base])
    for k in range(1,len(outside)+1):
        for T in itertools.combinations(outside,k):
            available+=(-1)**k*intersect_area([base,*[excluded[i] for i in T]])
    assert available>0,(points,masses,S,radius,available)
    return {'S':list(S),'m_S':str(mass),'test_radius':str(radius),
            'last_critical_below_mass_radius':str(last),
            'exact_feasible_center_area':str(available),
            'capacity_equals_m_S':True}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--certificates',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    start=time.perf_counter()
    raw=args.certificates.read_bytes();data=json.loads(raw)
    rows=data['certificates']
    assert len(rows)==data['P_equals_mass_exact_dual_certificates']==100
    assert [r['trial'] for r in rows]==list(range(100))
    results=[]
    for cert in rows:
        points=[tuple(map(F,p)) for p in cert['points']]
        masses=list(map(F,cert['masses']))
        N=len(points)
        assert N==5==len(masses) and len(set(points))==N and all(m>0 for m in masses)
        active=cert['positive_dual_support']
        seen=set();loads=[F(0)]*N;objective=F(0);checks=[]
        for item in active:
            S=tuple(item['S']);z=F(item['z']);c=F(item['c_S'])
            assert S and all(isinstance(i,int) and 0<=i<N for i in S)
            assert len(set(S))==len(S) and S not in seen and z>0
            seen.add(S)
            assert 0<c<=sum(masses[i] for i in S)
            for i in S:loads[i]+=z
            objective+=z*c
            # Equality with total mass and source loads forces saturation;
            # validate it directly and independently establish its geometry.
            assert c==sum(masses[i] for i in S)
            checks.append(capacity_saturation_check(points,masses,S))
        assert loads==list(map(F,cert['source_loads'])) and all(v<=1 for v in loads)
        assert objective==F(cert['dual_objective'])==F(cert['total_mass'])==sum(masses)
        results.append({'trial':cert['trial'],'active_capacities':checks,
                        'exact_source_loads':list(map(str,loads)),
                        'exact_objective':str(objective),'P_equals_total_mass':True})
    out={'certificates_verified':len(rows),
         'positive_capacity_coefficients_independently_verified':sum(len(r['active_capacities']) for r in results),
         'arithmetic':'Fraction and integer isqrt','LP_solver_used':False,
         'probe_coefficient_function_called':False,'random_search_rerun':False,
         'geometric_check':'independent inclusion-exclusion for exact center area in the last open radius layer below sqrt(m(S))/2',
         'certificate_file_sha256':hashlib.sha256(raw).hexdigest(),
         'audit_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'python':sys.version,'platform':platform.platform(),
         'verified_utc':datetime.now(timezone.utc).isoformat(),
         'elapsed_seconds':time.perf_counter()-start,
         'command':[sys.executable,*sys.argv],'certificates':results}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='certificates'},indent=2))

if __name__=='__main__':main()
