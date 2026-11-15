"""Independent exact R52 input-path and source-phase audit; no bbox/sweep replay."""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime,timezone
import hashlib,json,math,os,platform,sys,time

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(",",":")).encode()).hexdigest()
def ww(x):return [F(p[2])for p in x["pts"]]
def norm(x):return {**{k:str(F(x[k]))for k in ("Px","Py","lam")},"pts":[[str(F(v))for v in p]for p in x["pts"]]}

def main():
 start=time.perf_counter();assert os.path.ismount("/Volumes/data")
 topic=Path(__file__).resolve().parent.parent;raw=Path("/Volumes/data/CodexTemp/T2/planar-centered-square-weak11");geo=Path("/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-r001-geometry/planar-centered-square-weak11/computation/round052")
 bp=topic/"computation/inputs/r008_final_0611.json";mp=topic/"computation/inputs/r014_final_0935.json";vp=topic/"computation/results/round017/search_ledger.json";op=raw/"round043/cross-group-refine/candidate_0195.json";fp=raw/"round048/mass-geometry/cross/candidate_0010.json";binary=raw/"round017/binary-mass/manifest.json";cross=raw/"round048/mass-geometry/cross"
 base,model,original,final=map(read,(bp,mp,op,fp));vectorsource=read(vp)[0];manifest=read(geo/"slice/manifest.json");handoff=read(geo/"handoff.json");files=[bp,mp,vp,op,fp,binary,cross/"manifest.json",cross/"results.json",geo/"slice/manifest.json",geo/"handoff.json"]
 assert sha(bp)=="6c2e5d2d8fba5063bd45ecd4467f8c74cca5c0edbef958f90146a103029d5910" and sha(op)=="97c70896d07591e65e690be51b08a02b439f070440d3748aad12d0b2ae140dff" and sha(fp)=="b4e2f0f4d8132741747a4569e727431889dcb521840dd053b9a5f87527275e68"
 assert sha(binary)==vectorsource["manifest_sha256"]==manifest["source_hashes"]["mass_vectors"] and read(binary)["selected_vectors"]==vectorsource["selected_vectors"]
 for key,p in [("base",bp),("original",op),("cross_manifest",cross/"manifest.json"),("cross_results",cross/"results.json"),("script",geo/"prepare_slice.py")]:assert manifest["source_hashes"][key]==sha(p)
 # Reconstruct all eight archived vectors from independent R14 group/anchor data.
 n=len(base["pts"]);assert n==51 and len(original["pts"])==len(final["pts"])==204;origw=ww(model);vectors=[]
 for v in vectorsource["selected_vectors"]:
  a=v["aliases"][0];members=[a["quadrant"]*n+j for j in a["members"]];anchor=a["quadrant"]*n+a["anchor"];group=sum(origw[j]for j in members);cap=min(group,origw[anchor]);vec=[F(0)]*204
  for j in members:vec[j]+=a["sign"]*F(a["eps"])*cap*origw[j]/group
  vec[anchor]-=a["sign"]*F(a["eps"])*cap
  assert vec==list(map(F,v["vector"])) and sum(vec)==0;vectors.append(vec)
 def mass(c):return [w+sum(k*v[i]for k,v in zip(c,vectors))for i,w in enumerate(origw)]
 assert mass([F(-1,8)]*2+[F(1)]*2+[F(42,5)]*4)==ww(original)
 assert mass([F(-13,200)]*2+[F(1)]*2+[F(39,5)]*4)==ww(final)
 velocity=[F(3,50)*(vectors[0][i]+vectors[1][i])-F(3,5)*sum((vectors[j][i]for j in range(4,8)),F())for i in range(204)]
 assert velocity==[b-a for a,b in zip(ww(original),ww(final))]==list(map(F,read(geo/"slice/velocity_mass.json"))) and sum(velocity)==0
 assert all(a[:2]==b[:2]for a,b in zip(original["pts"],final["pts"]))
 P=[F(base[k])for k in ("Px","Py")];period=[F(original[k])for k in ("Px","Py")];lam=F(original["lam"]);total=sum(ww(original));z0=P[0]/1000000;zmax=P[0]/500000;H=P[0]/100000000
 assert period==[2*p for p in P] and lam==1 and total==sum(ww(final))==4*sum(ww(base))
 vx=read(geo/"slice/velocity_x.json");labels=[i for i,x in enumerate(vx)if x];assert labels==[130,131] and vx==[int(i in (130,131))for i in range(204)]
 assert labels==[2*n+j for j in (28,29)] and divmod(2,2)==(1,0)
 # Explicit unwrapped lift from R8: fixed shortest lift at t=0, then +z.
 deltas=[]
 for i,p in enumerate(original["pts"]):
  cx,cy=divmod(i//n,2)
  for axis,copy in ((0,cx),(1,cy)):
   d=(F(p[axis])-F(base["pts"][i%n][axis])-copy*P[axis])%period[axis]
   if d>period[axis]/2:d-=period[axis]
   deltas.extend([abs(d),abs(d+(zmax if axis==0 and vx[i] else 0))])
 delta=max(deltas);kappa=max(max(a,b)/F(base["pts"][i%n][2])for i,(a,b)in enumerate(zip(ww(original),ww(final))))
 assert kappa<=F(1001,1000) and min(ww(original)+ww(final))>0
 tail5=kappa*sum(ww(base))/(P[0]*P[1])*(1+(P[0]+2*delta)/10)*(1+(P[1]+2*delta)/10)
 tail10=total/(period[0]*period[1])*(1+period[0]/20)*(1+period[1]/20)
 assert tail5<lam and tail10<lam
 corners=[];ct={t["id"]:t for t in read(cross/"manifest.json")["tasks"]};cr={r["id"]:r for r in read(cross/"results.json")}
 for c in manifest["corner_certificates"]:
  t=F(c["t"]);zn=F(c["normalized_z"]);data=norm(original)
  for i,p in enumerate(data["pts"]):p[2]=str(F(p[2])+t*velocity[i]);p[0]=str((F(p[0])+vx[i]*zn*P[0])%period[0])
  dest=geo/f"slice/corner_t{int(t)}_z{int(zn!=0)}.json";assert norm(read(dest))==data==ct[c["task_id"]]["input"] and digest(data)==c["canonical_sha256"]==ct[c["task_id"]]["input_sha256"]==cr[c["task_id"]]["input_sha256"] and sha(dest)==c["file_sha256"]
  corners.append({"t":str(t),"z_physical":str(zn*P[0]),"task_id":c["task_id"],"canonical_sha256":digest(data),"file_sha256":sha(dest)});files.append(dest)
 mid=norm(original)
 for i,p in enumerate(mid["pts"]):p[2]=str(F(p[2])+velocity[i]/2)
 assert norm(read(geo/"slice/midmass_t_half.json"))==mid
 # This input has NOT already received z0; the slab applies x += vx*z0 once.
 assert all(F(p[0])==F(o[0])for p,o in zip(mid["pts"],original["pts"]))
 slices=[]
 for dirname,t in [("one_axis_t0",F(0)),("one_axis",F(1,2)),("one_axis_t1",F(1))]:
  sp=geo/dirname/"slab.json";s=read(sp);data=norm(original)
  for i,p in enumerate(data["pts"]):p[2]=str(F(p[2])+t*velocity[i])
  assert s["schema"]=="position-slab-enumerator-v1" and s["status"]=="complete" and s["periodic"] is True and s["input_canonical_sha256"]==digest(data) and s["direction_sha256"]==digest(vx) and s["velocities"]==vx and F(s["s0"])==z0 and F(s["requested_H"])==H and F(s["H"])>0
  assert s["domain"]==list(map(str,[F(0),period[0],F(0),period[1]])) and F(s["tail"]["bound"])==tail10 and F(s["tail"]["radius"])==10 and F(s["tail"]["mass"])==total and F(s["tail"]["lambda"])==lam
  assert digest(s["source_catalog"])==s["source_catalog_sha256"] and digest([s["x_classes"],s["y_classes"]])==s["source_order_sha256"]
  for p in s["source_catalog"]:
   i=p["label"];jx,jy=p["period_shift"];assert 0<=i<204 and p["v"]==vx[i] and F(p["w"])==F(data["pts"][i][2]) and F(p["x"])==F(data["pts"][i][0])+vx[i]*z0+jx*period[0] and F(p["y"])==F(data["pts"][i][1])+jy*period[1]
  slices.append({"mass_t":str(t),"input_canonical_sha256":digest(data),"source_catalog_count":len(s["source_catalog"]),"source_catalog_sha256":s["source_catalog_sha256"],"source_order_sha256":s["source_order_sha256"],"slab_file_sha256":sha(sp),"tail10":str(tail10),"source_geometry_phase":"unshifted input + vx*z0 + period_shift*target_period"});files.append(sp)
 assert len({s["source_order_sha256"]for s in slices})==1
 # Pin handoff existing files; large slabs were read once for the narrow phase audit.
 handoffchecks=[]
 for h in handoff["sources"]:
  p=geo/h["path"];assert p.stat().st_size==h["bytes"] and sha(p)==h["sha256"];handoffchecks.append(h);files.append(p)
 q5=math.ceil(F(5)/min(period))+1;q10=math.ceil(F(10)/min(period))+1
 out={"status":"passed","scope":"R52 exact R43-to-R48 mass path, q2 source identity/physical phase, current frozen handoff hashes and uniform tail hypotheses only; no bbox, endpoint-order, effective-wall, mixed-derivative or area replay","parameter_domain":{"t":["0","1"],"z_physical":["0",str(zmax)],"z_normalized_to_R8_Px":["0","1/500000"],"slice_center_z":str(z0),"slice_requested_H":str(H)},"R8_period_x":str(P[0]),"target_period_x":str(period[0]),"target_period_y":str(period[1]),"mass_total":str(total),"mass_derivative_sum":"0","mass_parameters_t0":{"a":"-1/8","b":"42/5"},"mass_parameters_t1":{"a":"-13/200","b":"39/5"},"mass_derivative_formula":"3/50*(V0+V1)-3/5*(V4+V5+V6+V7)","q2":{"quadrant_index_zero_based":2,"copy_x":1,"copy_y":0,"original_R8_labels_zero_based":[28,29],"target_labels_zero_based":labels,"physical_velocity":1},"corners":corners,"slices":slices,"uniform_domain_certificate":{"strict_positive_mass_min":str(min(ww(original)+ww(final))),"actual_mass_constant":True,"kappa_endpoint_max":str(kappa),"lifted_delta_bound":str(delta),"dominated_tail_R5":str(tail5),"dominated_tail_R5_margin":str(lam-tail5),"periodic_count_tail_R10":str(tail10),"periodic_count_tail_R10_margin":str(lam-tail10),"q_R5":q5,"q_R10":q10,"L1_inheritance":"Any separately certified inner area at any parameter in this domain transfers by (1-2q/N)^2 using actual mass; D=2^bits and epsilon=1/(2D). No new halo area is certified here."},"handoff_files_checked":handoffchecks,"files_sha256":{str(p):sha(p)for p in dict.fromkeys(files)},"script_sha256":sha(Path(__file__)),"python":sys.version,"platform":platform.platform(),"command":[sys.executable,*sys.argv],"verified_utc":datetime.now(timezone.utc).isoformat(),"elapsed_seconds":time.perf_counter()-start}
 dest=topic/"computation/results/round052/input-path-review.json";dest.parent.mkdir(parents=True,exist_ok=True);assert not dest.exists();dest.write_text(json.dumps(out,indent=2)+"\n");print(json.dumps({k:v for k,v in out.items()if k not in ("handoff_files_checked","files_sha256")},indent=2));print("output",dest,"sha256",sha(dest))
if __name__=="__main__":main()
