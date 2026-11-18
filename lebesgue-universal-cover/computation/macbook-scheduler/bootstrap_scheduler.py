"""Attach the known completed/in-flight campaigns and seed approved project work."""
import hashlib
import json
from pathlib import Path
import shutil
import sys

from scheduler import save

def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    root=Path(sys.argv[1]).resolve()
    base=root.parent
    project=base.parent/"lebesgue-tomar-v3-20261118"
    term=project/"terminal-083741-sol"
    outer=project/"outer-083741-sol"
    python=str(term/"pure-venv/bin/python")
    for name in ("inbox","runs","inputs","gpu-spool"):
        (root/name).mkdir(parents=True,exist_ok=True)
    if not (root/"config.json").exists():
        save(root/"config.json",{"cpu_budget":14,"metal_source":str(base/"seven.metal")})
    external=[]
    def attach(name,pid,command,cpu,result,gpu=False,exit_path=None):
        row={"id":name,"pid":pid,"command_contains":command,"cpu":cpu,"gpu":gpu,"result":result}
        if exit_path:
            row["exit_path"]=str(exit_path)
        external.append(row)
    seeds=json.loads((base/"seeds.json").read_bytes())
    for i in (0,1):
        path=base/f"queue-{i}/best.jsonl"
        attach(f"legacy-seven-gpu-{i}",86742+i,f"campaign.py gpu {i}",1,
               {"path":str(path),"type":"jsonl","rows":len(seeds[i::2]),
                "row_keys":["part","pos","score","epoch","gpu_seconds","pose"]},gpu=True)
    attach("legacy-seven-cpu",86744,"campaign.py cpu",15,
           {"path":str(base/"cpu-leaves.jsonl"),"type":"jsonl","rows":len(seeds),"row_keys":["part","pos","pose","score"]})
    attach("legacy-r165-outer",86345,"replay_tomar_083741_scheduled.py",15,
           {"path":str(outer/"full-work/stage-ledger.json"),"type":"json","equals":{"complete":True,"parts":1024}},
           exit_path=outer/"full.exit")
    attach("legacy-r165-ordinary",86691,"replay_finite_point_083741.py",15,
           {"path":str(term/"regular-work/replay.json"),"type":"json","equals":{"ordinary_complete":True,"terminal_trees":2049}},
           exit_path=term/"regular.exit")
    attach("legacy-r165-dynamic",86812,"pure",10,
           {"path":str(term/"pure.stdout.json"),"type":"json","equals":{"complete":True,"leaves":34}},exit_path=term/"pure.exit")
    save(root/"external.json",{"jobs":external})

    def lines(path):
        with path.open() as src:
            return [json.loads(s) for s in src]
    tasks=[]
    if all(not (base/f"queue-{j}/best.jsonl").stat().st_size==0 for j in (0,1)):
        by_lane=[lines(base/f"queue-{j}/best.jsonl") for j in (0,1)]
        gpu=[r for lane in by_lane for r in lane]
        cpu=lines(base/"cpu-leaves.jsonl")
        if all(len(v)==len(seeds[j::2]) for j,v in enumerate(by_lane)) and len(cpu)==len(seeds):
            verify_version="arb-support-2048/"+hash_file(base/"exact_pose.py")[:12]
            for family,rows in (("gpu",gpu),("cpu",cpu)):
                chosen={}
                for r in sorted(rows,key=lambda r:r["score"]):
                    chosen.setdefault(tuple(r["pose"]),r)
                group=root/"inputs"/f"strict-seven-{family}"
                group.mkdir(exist_ok=True)
                for f in ("exact_pose.py","finite-point-leaf-input.json","multi-host-points.json"):
                    shutil.copyfile(base/f,group/f)
                source=group/"frozen-candidates.jsonl"
                source.write_text("".join(json.dumps(r,separators=(",",":"))+"\n" for r in list(chosen.values())[:300]))
                for i in range(8):
                    tasks.append(dict(kind="cpu",stage="strict-seven-candidates",version=verify_version,
                        input_sha256=hash_file(source),scope=f"{family}: ranks {i}+8k below 300",priority=85,cpu=1,
                        argv=[python,"-I","-B",str(group/"exact_pose.py"),str(i)],cwd=str(group),
                        result={"path":"result.jsonl","type":"jsonl","stdout":True,"rows":len(range(i,min(300,len(chosen)),8)),
                                "row_keys":["index","pose","angle_count","outer_area"]}))
            data=(base/"inputs.txt").read_text().split()
            assert len(data)==96
            best=sorted(gpu,key=lambda r:r["score"])
            unique=[]
            used=set()
            for r in best:
                place=(r["part"],r["pos"])
                if place not in used:
                    unique.append(r)
                    used.add(place)
                if len(unique)>=12:
                    break
            native=root/"inputs"/"seven-incumbent-candidates"
            native.mkdir(exist_ok=True)
            for i,r in enumerate(unique):
                inp=native/f"start-{i:02}.txt"
                inp.write_text(" ".join([*map(str,r["pose"][:14]),*data[14:]])+"\n")
                tasks.append(dict(kind="gpu",driver="seven",stage="seven-local-incumbents",version="resident-seven-v1/"+hash_file(base/"seven.metal")[:12],
                    input=str(inp),input_sha256=hash_file(inp),scope=f"refine GPU winner from {r['part']}/{r['pos']}; independent seed {202611190+i}",
                    priority=25,cpu=1,epochs=12288,seed=202611190+i,
                    result={"path":"result.json","type":"json","equals":{"complete":True,"backend":"Metal"}}))
    contra=base/"counterexample-v1"
    if (contra/"run_remote.sh").is_file():
        checksum=hash_file(contra/"native-input.json")
        assert checksum=="17cd25f8c3fdda8989e99634174dccc1a4278701fcf3b7dc3b03b2ed1fcb31e5"
        tasks.append(dict(kind="gpu",stage="counterexample-scout",version="v1",
                          input_sha256=checksum,scope="19 anchors and 85 nearby hosts, both hands; shared-translation scoring",
                          priority=50,cpu=1,argv=["/bin/zsh",str(contra/"run_remote.sh"),"{run}/result.json","2048"],
                          cwd=str(contra),result={"path":"result.json","type":"json","equals":{"complete":True}}))
    assert tasks
    initial=root/"inbox"/"initial-20261118.json"
    if not initial.exists():
        save(initial,{"tasks":tasks})
    print(len(tasks), "initial tasks")

if __name__=="__main__":
    main()
