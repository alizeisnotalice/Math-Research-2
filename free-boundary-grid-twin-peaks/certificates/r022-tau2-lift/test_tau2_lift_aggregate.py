#!/usr/bin/env python3
"""Hand-weighted and adversarial tests for aggregate_tau2_lift.py."""
import hashlib,json,pathlib,subprocess,sys,tempfile
ROOT=pathlib.Path(__file__).parent; AGG=ROOT/'aggregate_tau2_lift.py'; BASE=0xA022000000000000
def sh(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def record(r,z,W,D,F,G,logz=None):
 n=3;P=2;seed=(BASE+n*1_000_000+20_000+r)&((1<<64)-1); means=[W/P,D/P,F/P,G/P]
 if logz is None: logz=__import__('math').log(z)
 zz=0 if logz < -744.4400719213812 else z
 return {'kind':'replicate','tau':2,'replicate':r,'seed_hex':f'{seed:016x}','Z':zz,'U':1,'logZ':logz,'U_over_Z':1,'D':1,'F':4,'G':0,'F_over_Z':36,'G_over_Z':0,'F_minus_2G_minus_4D':0,'dead':False,'tau2_lift_L_groups':[W],'tau2_lift_L_distance2_groups':[D],'tau2_lift_L_F0_groups':[F],'tau2_lift_L_G0_groups':[G],'tau2_lift_group_LF_minus_2LG_minus_4LD':[0],'tau2_lift_scaled_mean_L':means[0],'tau2_lift_scaled_mean_L_distance2':means[1],'tau2_lift_scaled_mean_L_F0':means[2],'tau2_lift_scaled_mean_L_G0':means[3],'tau2_lift_complete_K2_Z':zz*means[0],'tau2_lift_complete_K2_D':zz*means[1],'tau2_lift_complete_K2_F':zz*means[2],'tau2_lift_complete_K2_G':zz*means[3]}
def run(raw,side,ok):
 c=subprocess.run([sys.executable,str(AGG),str(raw),'--sidecar',str(side),'--source',str(AGG)],text=True,capture_output=True)
 if (c.returncode==0)!=ok:raise AssertionError((ok,c.stdout,c.stderr))
 return c
def main():
 meta={'kind':'metadata','implementation':'r022-tau2-terminal-lift-complete-K2','base_source_sha256':'x','n':3,'particles':2,'replicates':2,'tau_lo':2,'tau_hi':2,'domain':'0','method':'strat','move_steps':0,'lift_scope':'x','lift_integer_bound':'x'}
 end={'kind':'run_end','wall_seconds':1}; base=[meta,record(0,2,3,6,24,0),record(1,4,4,10,36,-2),end]
 with tempfile.TemporaryDirectory() as td:
  td=pathlib.Path(td)
  def put(name,rows,wrongmeta=False):
   raw=td/(name+'.jsonl');raw.write_text('\n'.join(json.dumps(x) for x in rows)+'\n'); side=td/(name+'.side.json');em={'n':3,'particles':2,'replicates':2,'tau_lo':2,'tau_hi':2,'domain':'0','method':'strat','move_steps':0};em['n']=4 if wrongmeta else em['n'];side.write_text(json.dumps({'raw_sha256':sh(raw),'source_sha256':sh(AGG),'source_path':str(AGG),'command':'synthetic','expected_metadata':em}));return raw,side
  raw,side=put('good',base);out=json.loads(run(raw,side,True).stdout);assert out['Q'].startswith('0.2626262626') and out['V'].startswith('1.090909') and out['C'].startswith('-0.090909')
  for name,mut in [('missing',lambda x:x.pop(2)),('seed',lambda x:x[2].update(seed_hex=x[1]['seed_hex'])),('false_dead',lambda x:x[1].update(dead=True)),('wronggroup',lambda x:x[1]['tau2_lift_L_groups'].__setitem__(0,1)),('fakeinteger',lambda x:x[1]['tau2_lift_L_groups'].__setitem__(0,'3.5'))]:
   rows=json.loads(json.dumps(base));mut(rows);r,s=put(name,rows);run(r,s,False)
  # This one is canonical-looking (Z=0, every target field null), so merely
  # checking payload consistency would wrongly accept it.
  rows=json.loads(json.dumps(base));x=rows[1];x['dead']=True;x['Z']=0
  for k in ('logZ','U_over_Z','D','F','G','F_over_Z','G_over_Z','F_minus_2G_minus_4D','tau2_lift_L_groups','tau2_lift_L_distance2_groups','tau2_lift_L_F0_groups','tau2_lift_L_G0_groups','tau2_lift_group_LF_minus_2LG_minus_4LD','tau2_lift_scaled_mean_L','tau2_lift_scaled_mean_L_distance2','tau2_lift_scaled_mean_L_F0','tau2_lift_scaled_mean_L_G0','tau2_lift_complete_K2_Z','tau2_lift_complete_K2_D','tau2_lift_complete_K2_F','tau2_lift_complete_K2_G'):x[k]=None
  r,s=put('canonical_dead',rows);run(r,s,False)
  rows=json.loads(json.dumps(base));rows[1]['tau2_lift_complete_K2_D']=0;r,s=put('falsezero',rows);run(r,s,False)
  r,s=put('wrongmetadata',json.loads(json.dumps(base)),True);run(r,s,False)
  rows=json.loads(json.dumps(base));rows[1].update(logZ=-745.0,Z=0,tau2_lift_complete_K2_Z=0,tau2_lift_complete_K2_D=0,tau2_lift_complete_K2_F=0,tau2_lift_complete_K2_G=0);r,s=put('nearunderflow',rows);run(r,s,False)
  # Nonzero subnormal-scale complete fields must not be silently zeroed.
  tiny=1e-310; rows=[meta,record(0,2,3,6,24,0),record(1,tiny,4,10,36,-2,__import__('math').log(tiny)),end];r,s=put('tiny_correct',rows);run(r,s,True)
  rows=json.loads(json.dumps(rows));rows[2].update(tau2_lift_complete_K2_Z=0,tau2_lift_complete_K2_D=0,tau2_lift_complete_K2_F=0,tau2_lift_complete_K2_G=0);r,s=put('tiny_falsezero',rows);run(r,s,False)
  # At the least Double subnormal, a product is quantized; the exact allowed
  # absolute discrepancy is half a subnormal, rather than an arbitrary floor.
  least=float.fromhex('0x0.0000000000001p-1022'); rows=[meta,record(0,2,3,6,24,0),record(1,least,4,10,36,-2,__import__('math').log(least)),end];r,s=put('least_subnormal',rows);run(r,s,True)
  rows=json.loads(json.dumps(base));rows[1].update(logZ=-1000,Z=0,tau2_lift_complete_K2_Z=0,tau2_lift_complete_K2_D=0,tau2_lift_complete_K2_F=0,tau2_lift_complete_K2_G=0);r,s=put('underflow',rows);run(r,s,True)
 print('tau2_lift_guard_tests=PASS: hand, missing, seed, dead, integer, metadata, underflow, tiny nonzero and least-subnormal quantization')
if __name__=='__main__':main()
