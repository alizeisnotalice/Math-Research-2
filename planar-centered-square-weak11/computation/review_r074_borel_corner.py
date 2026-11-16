"""Independent exact geometry replay, no producer imports. Python stdlib only."""
from fractions import Fraction as F
from pathlib import Path
import hashlib, json, os, platform, sys, time

start = time.monotonic()
source = Path(sys.argv[1]); outpath = Path(sys.argv[2])
raw = source.read_bytes(); claim = json.loads(raw)
n = claim['n']; c = F(claim['density']); N = claim['receiver_mesh']
e = F(1,n*n); h = F(n+1,n*n); w = 1-e; p = F(1,2)-e
coords = [i*h for i in range(n)]
def overlap(a,b,l,u): return max(F(0), min(b,u)-max(a,l))
def outside_area(x,y,r):
    return (overlap(x-r,x+r,-1,p)*overlap(y-r,y+r,-1,p)
            -overlap(x-r,x+r,-e,p)*overlap(y-r,y+r,-e,p))
def count(a,b): return sum(a <= s <= b for s in coords)
cuts = sorted({F(-1), w} | {z for a in coords+[-F(1),p,-e]
                                      for z in (a,a-1) if -1 <= z <= w})
assert cuts == list(map(F,claim['unit_origin_cuts']))
vertices = {(F(row['origin'][0]),F(row['origin'][1])):row
            for row in claim['unit_vertices']}
minslack = F(1); ratios = []
for u in cuts:
    for v in cuts:
        k = count(u,u+1)*count(v,v+1)
        area = outside_area(u+F(1,2),v+F(1,2),F(1,2))
        slack = 1-k*e-c*area
        row = vertices[u,v]
        assert (k,area,slack) == (row['internal_count'],F(row['outside_area']),F(row['slack']))
        assert slack >= 0
        minslack = min(minslack,slack)
        if area: ratios.append((1-k*e)/area)
assert len(vertices) == len(cuts)**2 == claim['unit_vertex_count']
assert min(ratios) == c
# All point/open-interval strata. On each open stratum the atom count is fixed;
# the continuous area is bilinear, so its supremum is at endpoint limits.
strata = [(x,x) for x in cuts] + list(zip(cuts,cuts[1:]))
limit_tests = 0
for a,b in strata:
    for d,f in strata:
        k = count((a+b)/2,(a+b)/2+1)*count((d+f)/2,(d+f)/2+1)
        for u in set((a,b)):
            for v in set((d,f)):
                assert k*e+c*outside_area(u+F(1,2),v+F(1,2),F(1,2)) <= 1
                limit_tests += 1
step = (F(1,2)-e)/N; seen=set(); minmargin=F(1); maxr=F(0)
for row in claim['cell_certificate']:
    i,j = row['i'],row['j']; assert 0 <= i < N and 0 <= j < N
    assert (i,j) not in seen; seen.add((i,j))
    x0 = -F(1,2)+i*step; x1=x0+step
    y0 = -F(1,2)+j*step; y1=y0+step; r=F(row['closed_radius'])
    assert 0 < r < F(1,2)
    # Literal atom distances to the entire receiver rectangle, no prefix rows.
    k = sum(max(abs(s-x0),abs(s-x1),abs(t-y0),abs(t-y1)) <= r
            for s in coords for t in coords)
    assert k == row['guaranteed_internal_count']
    assert x0-r > -1 and y0-r > -1 and x1+r < p and y1+r < p
    # Left/bottom stay below J and right/top inside J for every center in cell;
    # J-overlap is increasing in each center coordinate, minimum at lower corner.
    assert x1-r < -e and y1-r < -e
    assert x0+r >= -e and y0+r >= -e and x1+r <= 1 and y1+r <= 1
    minoutside = outside_area(x1,y1,r)
    margin = k*e+c*minoutside-4*r*r
    assert margin == F(row['strict_margin']) > 0
    minmargin=min(minmargin,margin); maxr=max(maxr,r)
area = len(seen)*step*step
assert len(seen) == claim['accepted_closed_cells']
assert area == F(claim['lower_area']) > F((n//2)**2,2*n*n)
excluded=set(); upper_tests=0
for row in claim.get('excluded_cell_certificate',[]):
    i,j=row['i'],row['j']; assert 0 <= i < N and 0 <= j < N
    assert (i,j) not in seen and (i,j) not in excluded; excluded.add((i,j))
    x0=-F(1,2)+i*step; x1=x0+step; y0=-F(1,2)+j*step; y1=y0+step
    bounds=[]
    for k in range(1,n+1):
        for l in range(1,n+1):
            # Any actual witness capturing k by l old sources contains this
            # prefix. Minimum required radius anywhere in the receiver cell.
            rmin=max(coords[k-1]-x1,coords[l-1]-y1)
            if rmin >= F(1,2): continue
            jmin=overlap(x0-rmin,x0+rmin,-e,1)*overlap(y0-rmin,y0+rmin,-e,1)
            bound=k*l*e-4*(1-c)*rmin*rmin-c*jmin
            assert bound <= 0
            bounds.append(bound); upper_tests += 1
    assert (str(max(bounds)) if bounds else None) == row['maximum_prefix_upper']
upper_area=(N*N-len(excluded))*step*step
assert upper_area == F(claim['upper_area']) < F(1,8)+F(3,4*n**4)
result = {'status':'PASS','method':'literal atom maxdist plus exact P/J overlap; all unit-origin point/open-interval strata',
          'n':n,'density':str(c),'unit_vertices':len(vertices),'unit_strata':len(strata)**2,
          'unit_stratum_corner_limits':limit_tests,'unit_minimum_slack':str(minslack),
          'accepted_cells':len(seen),'lower_area':str(area),'minimum_strict_margin':str(minmargin),
          'maximum_certificate_radius':str(maxr),'mass_half':str(F((n//2)**2,2*n*n)),
          'excluded_cells':len(excluded),'independent_upper_prefix_tests':upper_tests,
          'upper_area':str(upper_area),
          'lower_minus_mass_half':str(area-F((n//2)**2,2*n*n)),
          'scope':'Borel density only. Global unit completeness uses bilinear endpoint argument; unit-to-all-square uses density<=1. No finite atomic transfer and no main corner bound refutation.',
          'producer_result':str(source),'producer_result_sha256':hashlib.sha256(raw).hexdigest(),
          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'runtime':{'python':platform.python_version(),'pid':os.getpid(),'seconds':time.monotonic()-start}}
outpath.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result))
