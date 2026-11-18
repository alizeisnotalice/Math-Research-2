"""Keep a small seven-body search horizon queued around diverse genuine leaf starts."""
import hashlib
import json
from pathlib import Path
import sys

from scheduler import atomic, key_for, read_json, save

BATCH=8
EPOCHS=12288
SEED=202612000
CANDIDATE_LIMIT=2083

def digest(data):
    return hashlib.sha256(data).hexdigest()

def freeze(root):
    p=root/"inputs"/"seven-sweep.json"
    if p.exists():
        return p
    base=root.parent
    initial=read_json(root/"inbox"/"initial-20261118.json")["tasks"]
    already={g["scope"].split("from ",1)[1].split(";",1)[0]
             for g in initial if g.get("driver")=="seven"}
    buckets={}
    for lane in (0,1):
        for text in (base/f"queue-{lane}/best.jsonl").read_text().splitlines():
            r=json.loads(text)
            key=f"{r['part']}/{r['pos']}"
            if key not in already:
                buckets.setdefault(r["part"],[]).append(r)
    for seq in buckets.values():
        seq.sort(key=lambda r:(r["score"],r["pos"]),reverse=True)
    ordered=[]
    while any(buckets.values()):
        for part in sorted(buckets):
            if buckets[part]:
                ordered.append(buckets[part].pop())
    unique={f"{r['part']}/{r['pos']}":r for r in ordered}
    assert len(unique)==len(ordered)
    save(p,{"rows":ordered[:CANDIDATE_LIMIT]})
    return p

def manifest(root, index, wave):
    cutoff=wave*BATCH
    data=read_json(index)["rows"]
    root_dir=root/"inputs"/"seven-sweep"
    root_dir.mkdir(exist_ok=True)
    base=root.parent
    inputs=(base/"inputs.txt").read_text().split()
    assert len(inputs)==96
    version="seven-wide-refinement-v1/"+digest((base/"seven.metal").read_bytes())[:12]
    cpu_source=root/"inputs"/"strict-seven-gpu"
    verifier=root/"seven_followups.py"
    py=sys.executable
    jobs=[]
    parents=[]
    for i,r in enumerate(data[cutoff:cutoff+BATCH],start=cutoff):
        file=root_dir/f"row-{i:05}.txt"
        source=(" ".join([*map(str,r["pose"][:14]),*inputs[14:]])+"\n").encode()
        if file.exists():
            assert file.read_bytes()==source
        else:
            atomic(file,source)
        gpu=dict(kind="gpu",driver="seven",stage="seven-wide-refinement",version=version,input=str(file),
                 input_sha256=digest(source),scope=f"rank {i}, original {r['part']}/{r['pos']}, seed {SEED+i}",
                 priority=10,cpu=1,epochs=EPOCHS,seed=SEED+i,
                 result={"path":"result.json","type":"json","equals":{"complete":True,"backend":"Metal"}})
        jobs.append(gpu)
        parent=key_for(gpu)
        parents.append(parent)
        jobs.append(dict(kind="cpu",stage="strict-seven-Metal-followup",version="arb-wide-pose-v1",
            input_sha256=gpu["input_sha256"],scope=parent,priority=95,cpu=1,after=[parent],
            argv=[py,"-B",str(verifier),"verify",str(cpu_source),"{run}","native",parent],
            cwd=str(cpu_source),result={"path":"result.json","type":"json","equals":{"complete":True}}))
    if cutoff+BATCH<len(data):
        script=Path(__file__).resolve()
        jobs.append(dict(kind="cpu",stage="feed-seven-search",version="wave-v1/"+digest(script.read_bytes())[:12],
            scope=f"wave {wave+1}",priority=80,cpu=1,input_sha256=digest(index.read_bytes()),after=parents[:6],
            argv=[py,"-B",str(script),str(root),str(wave+1),"{run}/result.json"],cwd=str(root),
            result={"path":"result.json","type":"json","equals":{"complete":True}}))
    assert parents
    path=root/"inbox"/f"seven-sweep-{wave:04}.json"
    wanted={"tasks":jobs}
    if path.exists():
        assert read_json(path)==wanted
    else:
        save(path,wanted)
    return {"complete":True,"wave":wave,"native_batches":len(parents),"population":len(data),
            "manifest":path.name,"manifest_sha256":digest(path.read_bytes())}

if __name__=="__main__":
    root=Path(sys.argv[1]).resolve()
    reply=manifest(root,freeze(root),int(sys.argv[2]))
    if len(sys.argv)>3:
        save(Path(sys.argv[3]),reply)
    else:
        print(json.dumps(reply))
