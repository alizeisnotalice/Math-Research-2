#!/usr/bin/env python3
"""Exact checks for the cross-scale insertion counterexample, no search."""
from fractions import Fraction as F
from pathlib import Path
import argparse, json, platform, time

def frac(q): return str(q)
def certificate(n, epsilon):
    assert n >= 1 and 0 < epsilon <= 1
    points=[]; masses=[]; rows=[]; prefix=F(0); area=F(0)
    for i in range(n):
        p=F(2**i); mass=4*p*p-prefix; prefix+=mass
        assert mass > 0 and prefix == 4*p*p
        w=epsilon/(9*p); r=p+w
        slack=prefix+epsilon-4*r*r
        # Every center has -w<x<0 and p/2<abs(y)<p.
        # r is strictly larger than all distances to the origin and
        # earlier sources. Every later source is at least 2p away.
        assert slack > 0 and r < 2*p
        assert (p+w) == r and p < r and w < r
        if i: assert rows[-1]['p_exact'] == frac(p/2)
        cell_area=2*w*(p-p/2)
        assert cell_area == epsilon/9
        area+=cell_area
        points.append([frac(p),'0']); masses.append(frac(mass))
        rows.append(dict(i=i,p_exact=frac(p),prefix_mass=frac(prefix),
                         width=frac(w),radius=frac(r),strict_slack=frac(slack),
                         x_open=[frac(-w),'0'],
                         y_positive_open=[frac(p/2),frac(p)],
                         y_negative_open=[frac(-p),frac(-p/2)],
                         area=frac(cell_area)))
    assert prefix==4**n and area==n*epsilon/9
    eta=epsilon/(144*F(2**(n-1)))
    l1_area=F(0)
    for row in rows:
        p=F(row['p_exact']); w=F(row['width']); r=F(row['radius'])
        assert eta<=w/16 and eta<1/F(32)
        # Contracted center rectangle: -3w/4<x<-w/4,
        # 5p/8<abs(y)<7p/8. All mollified prefix mass fits.
        assert p+3*w/4+eta<r and 7*p/8+eta<r
        assert -w/4 < -eta
        # g's entire support stays too far vertically for density >1.
        assert 4*(5*p/8-eta)**2 > epsilon
        a=2*(w/2)*(p/4)
        assert a==epsilon/36
        l1_area+=a
    assert l1_area==n*epsilon/36
    return dict(n=n,epsilon=frac(epsilon),base_points=points,
                base_masses=masses,base_total=frac(prefix),
                inserted_point=['0','0'],inserted_mass=frac(epsilon),
                mollifier_square_halfside=frac(eta),l1_new_area_lower=frac(l1_area),
                disjoint_new_area_lower=frac(area),
                new_area_per_inserted_mass_lower=frac(area/epsilon),rows=rows)

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--out',required=True)
    args=parser.parse_args(); start=time.perf_counter()
    cases=[certificate(n,e) for n,e in [(1,F(1)),(2,F(1)),(8,F(1)),(32,F(1)),(32,F(1,1024))]]
    sources={'/Volumes/data/CodexProjects/T2/AGENTS.md': '78b0d90ef2438f81ccfc7a593cabdc3b4976deb7dbf9048b5c0426b82ff5e3a0', '/Volumes/data/CodexProjects/T2/PROJECTS.md': '687cc790f4438b0c66965cc21bb72e00ebc3502b2b2b1763159daf5378c4d582', '/Volumes/data/CodexProjects/T2/knowledge/README.md': '55dbe4de963c3b011f5ae05763d46218d2d29c4ddf84e25595ded8f3bb63717b', '/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/STATE.md': 'b1b68fb887782ebcf42e8684568fb7590b5d137cdb313dfa092e98730ad5f8bd', '/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/proofs/general-dual-interface.md': '8b9184ac14d3e38cb31c82e5c44a697e366a2c46cc3294323c8ffa393627c821', '/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/proofs/kernel-family-floor.md': '9c8d1f05b4578e77a94bcc9c6d0c8aff0a64bfa26e0621f1a9ad9537db93be15', '/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/rounds/015.md': '962de27e8b5a5173c964b7d995776284f196ce64fff007122d433b7724da7d1a', '/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-general-dual/planar-centered-square-weak11/proofs/mass-capacity-equality-attempt.md': '7cc74c83702818078187d9f1b29dc070e4473cbd7e0294fbc3b7934bde681740'}
    result=dict(status='passed',method='Python standard-library Fraction; explicit construction, no search',
      python=platform.python_version(),random_seed=None,tolerance=None,
      elapsed_seconds=time.perf_counter()-start,cases=cases,source_sha256=sources,
      limitations='Checks finite rational witness rectangles. All-N theorem and absence of baseline halo in x<=0 are proved analytically in the accompanying note; no complete halo area computed.')
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','python','elapsed_seconds']}))
    print('cases='+str(len(cases))+'; strongest explicit new-area/mass lower='+cases[3]['new_area_per_inserted_mass_lower'])
if __name__=='__main__': main()
