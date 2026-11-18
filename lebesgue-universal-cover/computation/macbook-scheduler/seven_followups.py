"""Add exact checks after actual Metal completions and independent backup shards."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

from scheduler import atomic, key_for, read_json, save

def sha(data):
    return hashlib.sha256(data).hexdigest()

def rows(path):
    with path.open() as f:
        return [json.loads(line) for line in f]

def verifier(source, run, kind, argv):
    sys.path.insert(0,str(source))
    from exact_pose import calculate
    from flint import arb
    if kind=="native":
        parent=argv[0]
        root=source.parents[1]
        meta=read_json(root/"ledger.json")["tasks"][parent]
        assert meta["state"]=="complete" and meta["task"]["driver"]=="seven"
        result=meta["result"]
        raw=Path(result["path"]).read_bytes()
        assert sha(raw)==result["sha256"]
        candidate=json.loads(raw)
        assert candidate["backend"]=="Metal" and candidate["complete"] and len(candidate["pose"])==18
        area=calculate(candidate,4096)
        assert arb(area)>0
        save(run/"result.json",{"complete":True,"mode":"fixed-seventh-Metal-pose",
            "parent":parent,"candidate_sha256":result["sha256"],"angle_count":4096,
            "outer_area":area,"pose":candidate["pose"]})
    else:
        index=Path(argv[0])
        first,last=map(int,argv[1:3])
        cloud=rows(index)
        assert 0<=first<last<=len(cloud)
        out=[]
        for k in range(first,last):
            candidate=cloud[k]
            area=calculate(candidate,1024)
            assert arb(area)>0
            out.append(dict(index=k,angle_count=1024,part=candidate["part"],pos=candidate["pos"],
                            outer_area=area,pose=candidate["pose"]))
        encoded="".join(json.dumps(x,separators=(",",":"))+"\n" for x in out)
        atomic(run/"result.jsonl",encoded.encode())

def producer(root):
    base=root.parent
    py=sys.executable
    source=root/"inputs"/"strict-seven-gpu"
    script=Path(__file__).resolve()
    policy="arb-separated-pose-certs/"+sha(script.read_bytes())[:12]
    init=read_json(root/"inbox"/"initial-20261118.json")["tasks"]
    tasks=[]
    for g in init:
        if g.get("driver")!="seven":
            continue
        parent=key_for(g)
        tasks.append(dict(kind="cpu",stage="strict-seven-Metal-followup",version=policy,
            input_sha256=g["input_sha256"],scope=parent,priority=95,cpu=1,after=[parent],
            argv=[py,"-B",str(script),"verify",str(source),"{run}","native",parent],
            cwd=str(source),result={"path":"result.json","type":"json","equals":{"complete":True}}))
    seen=set()
    families=[
        ("gpu",rows(base/"queue-0/best.jsonl")+rows(base/"queue-1/best.jsonl")),
        ("cpu",rows(base/"cpu-leaves.jsonl")),
    ]
    for name,_ in families:
        seen.update(tuple(r["pose"]) for r in rows(root/"inputs"/f"strict-seven-{name}"/"frozen-candidates.jsonl"))
    union=[]
    for name,seq in families:
        for r in sorted(seq,key=lambda r:r["score"]):
            pose=tuple(r["pose"])
            if pose not in seen:
                seen.add(pose)
                union.append(r)
    source_file=root/"inputs"/"backup-seven-distinct.jsonl"
    data="".join(json.dumps(x,separators=(",",":"))+"\n" for x in union).encode()
    atomic(source_file,data)
    for a in range(0,len(union),48):
        b=min(a+48,len(union))
        tasks.append(dict(kind="cpu",stage="strict-seven-backup",version=policy,input_sha256=sha(data),
            scope=f"remaining independent poses [{a},{b})",priority=5,cpu=1,
            argv=[py,"-B",str(script),"verify",str(source),"{run}","backup",str(source_file),str(a),str(b)],
            cwd=str(source),result={"path":"result.jsonl","type":"jsonl","rows":b-a,
                                     "row_keys":["index","angle_count","outer_area","pose"]}))
    save(root/"inbox"/"strict-continuation-20261118.json",{"tasks":tasks})
    print("dependent strict tasks",len(init)-17,"new unique backup poses",len(union),
          "total tasks",len(tasks))

def main():
    if sys.argv[1]=="produce":
        producer(Path(sys.argv[2]).resolve())
    else:
        assert sys.argv[1]=="verify"
        verifier(Path(sys.argv[2]),Path(sys.argv[3]),sys.argv[4],sys.argv[5:])

if __name__=="__main__":
    main()
