"""Rebuild the historical selection and audit or replay every fresh original-area obligation."""
import argparse,hashlib,json,subprocess,sys,tempfile
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from zipfile import ZipFile
from flint import arb,ctx
p=argparse.ArgumentParser(description=__doc__)
for name in ("upstream","certdir","history_archive","history_results","masks","mask_index","table","archive","scratch"):p.add_argument(name,type=Path)
p.add_argument("--replay-weak",action="store_true");a=p.parse_args()
here=Path(__file__).resolve().parent;sys.path.insert(0,str(here))
from verify_tomar_fan_table import verified_table
def sha_bytes(x):return hashlib.sha256(x).hexdigest()
def sha(p):return sha_bytes(p.read_bytes())
with tempfile.TemporaryDirectory(dir=a.scratch) as tmp, ZipFile(a.archive) as z:
 d=Path(tmp);raw=z.read("selection.json");selection=json.loads(raw);target=selection["target"]
 file=d/"selection.json";file.write_bytes(raw)
 args=[a.upstream,a.certdir,a.history_archive,a.history_results,a.masks,a.mask_index,d]
 subprocess.run([sys.executable,"-I","-B",here/"collect_tomar_historical_weak.py",*args,"--target",target,"--recompute-slice"],stdout=subprocess.DEVNULL,check=True)
 assert (d/("weak-"+target.replace(".","")+".json")).read_bytes()==raw
 assert selection["canonical_outer_leaves"]==6710314
 dense="fe5b816060d713642dc44e3f8981f4ae5d1b2b1c1e271b7dff75f031dd75d87c"
 assert sha(a.table)==dense
 meta=json.loads((a.certdir/"part0000.json").read_bytes());verified_table(a.table,a.certdir,meta)
 expected_files={"selection.json","deeper-selection.json"}
 deep_raw=z.read("deeper-selection.json");deep=json.loads(deep_raw);deep_file=d/"deeper-selection.json";deep_file.write_bytes(deep_raw)
 assert deep["target"]==target and deep["primary_selection_sha256"]==sha_bytes(raw)
 weak_rows=selection["weak"];lookup={(r[0],r[1]):r for r in weak_rows};assert len(lookup)==len(weak_rows)==selection["canonical_weak"]
 expected_deeper=[];parents={};levels=[]
 for level,filename,rec,limit,depth in (("first",file,selection,32767,45),("deeper",deep_file,deep,1048575,60)):
  selected=defaultdict(list)
  for r in rec["weak"]:selected[r[0]].append(r[:3])
  checked=closed=calls=0;worst=0;dst=d/("replay-"+level)
  if a.replay_weak:
   subprocess.run([sys.executable,"-I","-B",here/"check_tomar_selected_fans.py",a.certdir,a.table,here/"tomar_parabolic_fan.cpp",dst,target,
    "--selected",filename,"--calls",str(limit),"--depth",str(depth)],stdout=subprocess.DEVNULL,check=True)
  for part in sorted(selected):
   name=f"{level}/part{part:04}.json";expected_files.add(name);b=z.read(name);res=json.loads(b)
   assert res["part"]==part and res["target"]==target and res["selection_sha256"]==sha(filename)
   assert res["depth"]==depth and res["calls"]==limit and res["source_sha256"]==sha(here/"tomar_parabolic_fan.cpp")
   assert res["driver_sha256"]==sha(here/"check_tomar_selected_fans.py") and res["decoder_sha256"]==sha(here/"tomar_leaf_box.py")
   assert res["table_sha256"]==dense and res["library_sha256"]=="15a3f277a1efd69564b41ee89e80a1ea1aa37fdb980bce51146099e61858b57e"
   assert res["header_sha256"]==sha(a.certdir/f"part{part:04}.json") and res["codes_sha256"]==sha(a.certdir/f"part{part:04}.npy")
   assert [[part,*row[:2]] for row in res["leaves"]]==selected[part]
   assert len(res["leaves"])==res["checked"] and sum(x[-1] for x in res["leaves"])==res["closed"]
   if a.replay_weak:
    replay=json.loads((dst/f"part{part:04}.json").read_bytes());assert replay["leaves"]==res["leaves"] and replay["library_sha256"]==res["library_sha256"]
   if level=="first":parents[str(part)]=sha_bytes(b)
   for idx,tag,num,k,term,ok in res["leaves"]:
    assert ok in (0,1) and 1<=num<=limit and 0<=k<=depth and (not ok or term>=1)
    if level=="first" and not ok:expected_deeper.append(lookup[(part,idx)])
    if level=="deeper":assert ok==1
    calls+=num;worst=max(worst,k)
   checked+=res["checked"];closed+=res["closed"]
  levels.append({"level":level,"checked":checked,"closed":closed,"integer_box_calls":calls,"maximum_depth":worst})
 assert deep["weak"]==expected_deeper and deep["parents"]==parents
 assert len(z.namelist())==len(expected_files) and set(z.namelist())==expected_files
 ctx.prec=128;outside=[]
 for i,n in enumerate(meta["ns"]):
  j=3*i;r=1-1/(2*(arb.pi()/(2*n)).cos());dist=min(-meta["lo"][j],meta["hi"][j],-meta["lo"][j+1],meta["hi"][j+1])
  area=arb.pi()*(arb(1)/4+r*r)/2+arb(dist)*(arb(1)/2+r);assert area>arb(target);outside.append(str(area))
 assert sum(x["closed"] for x in levels)==selection["canonical_weak"]
 print(json.dumps({"target":target,"result_kind":"conditional on author historical strong-leaf area arithmetic",
   "historical_strong":selection["canonical_strong"],"fresh_weak":selection["canonical_weak"],"canonical_leaves":selection["canonical_outer_leaves"],
   "full_original_leaves":selection["full_leaves"],"selection_data":selection["data"],"levels":levels,"fresh_weak_replayed":a.replay_weak,
   "root_exterior":outside,"slice_index_sha256":selection["slice_index_sha256"],"archive_bytes":a.archive.stat().st_size,
   "archive_sha256":sha(a.archive)},indent=2))
