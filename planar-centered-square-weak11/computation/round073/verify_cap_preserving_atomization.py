from fractions import Fraction as F
from pathlib import Path
import json,sys,platform,time,os,hashlib

def main():
 start=time.time();counts={'positive_count_pairs':0,'radius_identity_cases':0,'zero_capture_cases':0}
 for n in [2,3,5,7,9,17,33]:
  t=F(n,n+1)
  for delta in [F(1,2*n),F(1,10*n)]:
   c=t/(1+2*delta)
   for k in range(1,n+1):
    a=F(k*(n+1)-n,n*n)
    for l in range(1,n+1):
     b=F(l*(n+1)-n,n*n);G=F((n-k)*(n-l),n**3)
     assert 1-t-(F(k*l,n*n)-t*a*b)==G
     for R in [F(1),F(3,2),F(3)]:
      excess=F(k*l,n*n)+c*((R+2*delta)**2-(a+2*delta)*(b+2*delta))-R*R
      formula=-(1-2*n*delta)*G/(1+2*delta)-(R-1)*((1-c)*(R+1)-4*c*delta)
      assert excess==formula and excess<=0
      assert (1-c)*(R+1)-4*c*delta==(1-c)*(R-1)+2*(1-t)
      counts['radius_identity_cases']+=1
     counts['positive_count_pairs']+=1
   for R in [F(1),F(3,2),F(3)]:
    assert c*(R+2*delta)**2<=R*R;counts['zero_capture_cases']+=1
 n=9;delta=F(1,100);lam=F(9999,10000);t=F(n,n+1);c=t/(1+2*delta);e=F(1,n*n);w=1-e
 floor=lambda z:z.numerator//z.denominator
 ceil=lambda z:-((-z.numerator)//z.denominator)
 def inside_axis(a,b):return (ceil(a/delta),floor(b/delta)-1)
 def intersects_axis(a,b):return (ceil(a/delta)-1,floor(b/delta))
 # The exact lower endpoint is ceil(a/delta-1), equal to ceil(a/delta)-1.
 def count(pair):return max(0,pair[1]-pair[0]+1)
 B=inside_axis(F(-1),F(2));C=inside_axis(F(-1,2),w+F(1,2));Jhit=intersects_axis(-e-delta,1+delta)
 assert B==(-100,199) and C==(-50,147) and Jhit==(-3,101)
 N=count(B)**2-count(Jhit)**2;Nring=count(C)**2-count(Jhit)**2
 atom_mass=lam*c*delta**2;total=1+N*atom_mass;span=F(299,100)
 lower=(1+e)**2+Nring*atom_mass
 assert lower>F(7,2)
 assert span**2-F(35,33)*total>0
 assert 2*atom_mass>delta**2
 assert F(2,n*n)>F(n+1,n*n)**2
 assert F(1,n*n)+atom_mass>F(7,200)**2
 out={'schema':'r073-exact-cap-repair-scalar-v1','counts':counts,'fixture':{'n':n,'delta':str(delta),'contraction':str(lam),'density_before_contraction':str(c),'B_cell_indices':B,'C_cell_indices':C,'closed_expanded_J_intersection_indices':Jhit,'outside_atom_count_formula':N,'total_source_count_formula':N+n*n,'ring_singleton_count_formula':Nring,'each_outside_atom_mass':str(atom_mass),'actual_total_mass':str(total),'actual_span':str(span),'nonterminal_gap':str(span**2-F(35,33)*total),'local_full_halo_area_lower':str(lower),'local_full_halo_area_lower_decimal':float(lower),'strict_excess_over_7_over_2':str(lower-F(7,2)),'outside_adjacent_pair_margin':str(2*atom_mass-delta**2),'internal_adjacent_pair_margin':str(F(2,n*n)-F(n+1,n*n)**2),'one_bridge_pair_margin':str(F(1,n*n)+atom_mass-F(7,200)**2)},'scope':'Exact scalar identity samples and formula-based rectangular cell counts; not enumeration of79056sources, all their subsets, or complete full/residual/deletion halos. General all-R/all-caps proof is separate. Strict lower area concerns Hfull intersect C, not actual deletion.','runtime':{'python':platform.python_version(),'pid':os.getpid(),'elapsed_seconds':time.time()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}}
 Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':'PASS','counts':counts,'fixture':out['fixture']}))
if __name__=='__main__':main()
