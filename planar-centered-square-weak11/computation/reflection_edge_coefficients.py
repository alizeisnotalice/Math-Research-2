#!/usr/bin/env python3
"""Exact interval ledger for bulk/edge/corner polynomials and true-axis exclusion."""
import argparse,hashlib,json,os,sys
from fractions import Fraction as F
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def linear(terms):
    return (sum(c*(z[0] if c>=0 else z[1]) for c,z in terms),sum(c*(z[1] if c>=0 else z[0]) for c,z in terms))
def rec(z):return {'lower':str(z[0]),'upper':str(z[1]),'lower_float':float(z[0]),'upper_float':float(z[1])}
def inflate(z,e):return z[0]-e,z[1]+e

def main():
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('--tag',required=True);ap.add_argument('--state-dir',required=True);ap.add_argument('--bulk',required=True);ap.add_argument('--axis-error',required=True);ap.add_argument('--output',required=True);args=ap.parse_args()
    assert os.path.ismount('/Volumes/data');src=Path(args.input);ih=sha(src);data=json.loads(src.read_text());assert F(data['lam'])==1;m=sum(F(p[2]) for p in data['pts']);states={};sources=[];directory=Path(args.state_dir)
    names=['left','right','bottom','top','seam-old','seam-new','corner-lb','corner-lt','corner-rb','corner-rt','seam-bottom-old','seam-bottom-new','seam-top-old','seam-top-new']
    for state in ['bulk']+names:
        p=Path(args.bulk) if state=='bulk' else directory/f'{args.tag}-{state}.json';r=json.loads(p.read_text())
        assert r['input_sha256']==ih
        if state!='bulk':assert r['status']=='PASS' and r['state']==state and F(r['tail_bound'])<F(data['lam'])
        z=F(r['inner_area']),F(r['outer_area']);assert z[0]<=z[1];states[state]=z
        sources.append({'state':state,'path':str(p),'sha256':sha(p),'size':p.stat().st_size,'area':rec(z),'script_sha256':r.get('script_sha256',r.get('code_sha256')),'counts':r['counts']})
    errpath=Path(args.axis_error);err=json.loads(errpath.read_text());assert err['status']=='PASS' and err['input_sha256']==ih
    eg,eG=F(err['gamma_change_upper']),F(err['Gamma_change_upper'])
    A=states['bulk'];E=linear([(1,states[s]) for s in ['left','right','bottom','top']]);C=linear([(1,states[s]) for s in ['corner-lb','corner-lt','corner-rb','corner-rt']]);g=linear([(1,states['seam-new']),(-1,states['seam-old'])]);G=linear([(1,states['seam-bottom-new']),(-1,states['seam-bottom-old']),(1,states['seam-top-new']),(-1,states['seam-top-old'])]);b=linear([(1,E),(-8,A)]);c=linear([(16,A),(-4,E),(1,C)]);L=linear([(1,E),(-8,A),(1,g)]);K=linear([(16,A),(-4,E),(1,C),(-6,g),(1,G)]);F6=linear([(2,E),(-32,A),(1,C),(1,G)])
    trueL=inflate(L,eg);trueK=inflate(K,6*eg+eG);trueF6=inflate(F6,eG)
    assert L[1]<0 and F6[1]<0 and trueL[1]<-9 and trueF6[1]<-53
    result={'status':'PASS','python':sys.version,'input_sha256':ih,'script_sha256':sha(__file__),'tag':args.tag,'mass_per_cell':str(m),'sources':sources,'axis_error_source':{'path':str(errpath),'sha256':sha(errpath),'gamma_error':str(eg),'Gamma_error':str(eG)},'coefficient_intervals':{name:rec(z) for name,z in [('A',A),('E',E),('C_corner',C),('gamma_rational_axis',g),('Gamma_rational_axis',G),('old_linear_b',b),('old_constant_c',c),('reflection_leading_L_rational',L),('reflection_constant_K_rational',K),('F6_rational',F6),('reflection_leading_L_true_axis',trueL),('reflection_constant_K_true_axis',trueK),('F6_true_axis',trueF6)]},'formulas':{'old_area_N_ge_4':'A_N=(N-4)^2*A+(N-4)*E+C_corner','added_area_rational_axis_N_ge_6':'D_N=2*(N-6)*gamma+2*Gamma','reflected_ratio_difference_N_ge_6':'J_ref-A/m=[L*(N-6)+F6]/(m*N^2)','true_axis':'Same polynomial structure; gamma and Gamma enclosed by the certified perturbation bounds.'},'all_integer_N_ge_6_strictly_below_exact_periodic_bulk':True,'clean_rational_bounds':{'L_true_upper_strict':'-9','F6_true_upper_strict':'-53','ratio_deficit_lower':'(9*N-1)/(m*N^2), every integer N>=6'},'uncovered_scope':['N=1,2,3,4,5 may have different local states and are not excluded here.','This is only the stated N-by-N canonical cell truncation and its rightmost reflection, not all finiteizations or all input families.','No new general upper bound on C2.'],'conclusion':'Both L_upper<0 and F6_upper<0 hold after transfer to the true irrational axis, so every integer N>=6 in the specified reflected finite family has weak ratio strictly smaller than this template exact infinite-periodic ratio.'}
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','tag':args.tag,'L_true':rec(trueL),'F6_true':rec(trueF6),'all_N_ge_6_excluded':True},ensure_ascii=False))
if __name__=='__main__':main()
