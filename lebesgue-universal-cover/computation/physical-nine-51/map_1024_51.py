"""Map all 1024 active shards to native and original-leaf .83751 receipts."""
import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
for name in ("tree","seed","exterior_audit","old_core","theorem","output"):
    p.add_argument(name,type=Path)
p.add_argument("--replay",type=Path,action="append",required=True)
p.add_argument("--direct-audit",type=Path,action="append",required=True)
a = p.parse_args()
read=lambda x:json.loads(x.read_bytes())
sha=lambda x:sha256(x.read_bytes()).hexdigest()
target=Fraction(".83751")
tree_manifest=a.tree/"SHA256SUMS"
manifest={n.removeprefix("./"):h for h,n in (l.split() for l in tree_manifest.read_text().splitlines())}
seed=read(a.seed)
assert seed["tree_manifest_sha256"]==sha(tree_manifest)
assert Fraction(seed["target"])==Fraction(seed["radial_target"])==target
assert seed["target_integer"]==int(target*2**57)+1
rows={key: (n,selected,index_sha,score_sha) for key,n,_,_,selected,index_sha,score_sha in seed["rows"]}
assert len(rows)==1024 and sum(x[0] for x in rows.values())==6710314
ext=read(a.exterior_audit)
assert ext["all_four_layers"] and ext["parts"]==1024 and len(ext["replays"])==len(a.replay)
final=read(a.theorem)
assert final["complete_global"] and Fraction(final["strict_target"])==target
prior={int(p):sorted(pos) for p,pos in read(a.old_core).items()}
sources_for_originals=defaultdict(dict)
for k,v in final["known_originals"].items():
    part,pos=map(int,k.split("/"))
    sources_for_originals[part][pos]=v
assert sum(map(len,sources_for_originals.values()))==5442
assert all(set(sources_for_originals[i])==set(prior[i]) for i in prior)
direct=defaultdict(dict)
for f in a.direct_audit:
    d=read(f)
    digest=sha(f)
    assert Fraction(d["target"])==target and not d["pending"] and d["closed"]==d["original_roots"]
    for k,r in d["per_root"].items():
        part,pos=map(int,k.split("/"))
        assert pos not in direct[part] and pos not in sources_for_originals[part]
        direct[part][pos]=dict(families=r["families"],visits=r["nodes"],audit_sha256=digest)
assert sum(map(len,direct.values()))==512
covered={}
for i,work in enumerate(a.replay):
    log=work/"stage-ledger.json"
    assert sha(log)==ext["replays"][str(i)]["ledger_sha256"]
    ledger=read(log)
    assert ledger["seed_receipt_sha"]==sha(a.seed)
    assert len(ledger["audits"])==4
    scopes=[dict(layer["work_files"]) for layer in ledger["audits"]]
    scope=set(scopes[0])
    assert all(set(s)==scope for s in scopes)
    assert tuple(sorted(scope))==tuple(ext["replays"][str(i)]["parts"])
    assert not scope&set(covered)
    retain=read(work/"retained.json")
    for part in sorted(scope):
        n,selected,hp,hf=rows[part]
        before=selected
        name=f"part{part:04}"
        logs=[]
        for k in range(1,5):
            path=work/f"stage{k}"/(name+".json")
            rec=read(path)
            assert rec["config"]["part"]==part and rec["config"]["layer"]==k
            assert rec["config"]["goal"]==seed["target_integer"] and Fraction(rec["config"]["radial_goal"])==target
            assert rec["data_sha256"]==scopes[k-1][part]
            assert rec["tested"]==before==rec["remaining"]+rec["accepted"]
            before=rec["remaining"]
            logs.append(dict(inputs=rec["tested"],accepted_whole_leaves=rec["accepted"],remaining=before,
                             calls=rec["counts"][0],native_npz_sha256=rec["data_sha256"],
                             native_data_receipt_sha256=sha(path)))
        actual=sorted(set(direct.get(part,{}))|set(prior.get(part,[])))
        assert actual==retain.get(str(part),[]) and before==len(actual)
        assert n==n-selected+sum(r["accepted_whole_leaves"] for r in logs)+len(actual)
        for suffix in (".json",".npy"):
            f=a.tree/(name+suffix)
            assert sha(f)==manifest[name+suffix]
        covered[part]=dict(canonical_original_leaves=n,seed_count=selected,
            before_new_native_screen=n-selected,seed_original_positions_sha256=hp,seed_lower_scores_sha256=hf,
            geometry_component=i,original_header_sha256=manifest[name+".json"],
            original_preorder_array_sha256=manifest[name+".npy"],native_stages=logs,
            reopened_exterior=direct.get(part,{}),old_core_originals=sources_for_originals.get(part,{}))
assert set(covered)==set(rows)
result=dict(target=str(target),area_grid_denominator=2**57,integer_goal=seed["target_integer"],
            sha256_tree_manifest=sha(tree_manifest),sha256_theorem=sha(a.theorem),sha256_seed=sha(a.seed),
            canonical_shards=len(covered),canonical_original_leaves=sum(r["canonical_original_leaves"] for r in covered.values()),
            old_core_originals=sum(len(r["old_core_originals"]) for r in covered.values()),
            reopened_exterior_originals=sum(len(r["reopened_exterior"]) for r in covered.values()),
            shards={i:r for i,r in sorted(covered.items())})
a.output.write_text(json.dumps(result,separators=(",",":"))+"\n")
print(json.dumps({k:result[k] for k in ("canonical_shards","canonical_original_leaves","old_core_originals","reopened_exterior_originals")}))
