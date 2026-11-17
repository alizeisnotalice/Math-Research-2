#!/usr/bin/env python3
"""Build a compact, byte-bound archive index for round-055 degree-8 S={1}.

This is deliberately an archival transformer, not a mathematical re-run.  It
reads frozen JSON outputs, writes canonical JSONL/gzip material, and records
both archive hashes and the source-byte hashes from which it was made.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path

SELF = Path(__file__).resolve().parent
ENTRY = Path("/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/section_eigen39/round055-degree8-entry")
ALG = Path("/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/algebraic_global/round055-s1-independent")
OUT = SELF / "package"

def sha_path(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()

def compact(o) -> bytes:
    return json.dumps(o, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode() + b"\n"

def put_json(rel: str, obj) -> Path:
    p = OUT / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(compact(obj))
    return p

def put_gz_jsonl(rel: str, rows: list[object]) -> tuple[Path, str]:
    """Write deterministic gzip (mtime=0); return byte path and raw SHA."""
    p = OUT / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    raw_hash = hashlib.sha256()
    with p.open("wb") as fp:
        with gzip.GzipFile(filename="", mode="wb", fileobj=fp, mtime=0) as gz:
            for row in rows:
                line = compact(row)
                gz.write(line)
                raw_hash.update(line)
    if p.stat().st_size >= 10_000_000:
        raise RuntimeError(f"archive part too large: {p} ({p.stat().st_size})")
    return p, raw_hash.hexdigest()

def copy_frozen(source: Path, rel: str) -> dict:
    if source.is_symlink():
        raise RuntimeError(f"refuse symlink source {source}")
    target = OUT / rel
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target, follow_symlinks=False)
    if target.is_symlink() or sha_path(source) != sha_path(target):
        raise RuntimeError(f"copy mismatch {source} -> {target}")
    return {"archive": rel, "source": str(source), "sha256": sha_path(source), "bytes": source.stat().st_size}

def source_ledger(directory: Path) -> list[dict]:
    rows=[]
    for p in sorted((directory / "cores").glob("*.json")):
        x=json.loads(p.read_text())
        item={"name":p.name,"sha256":sha_path(p),"bytes":p.stat().st_size,
              "core":x["binding"]["core"],"status":x["status"]}
        if "records" in x: item["records"]=len(x["records"])
        else:
            item["row4_records"]=len(x["row4_records"])
            item["row7_records"]=len(x["row7_records"])
        rows.append(item)
    return rows

def h_from_row(row: list[int]) -> list[int]:
    T,N,b,c,d,e,f,g = row
    F=[1,0,b,c,d,e,f,g,1]
    h=[0]*17
    for i,x in enumerate(F):
        for j,y in enumerate(F): h[i+j]+=x*y
    for i,x in enumerate(F): h[i+1]+=T*x
    h[2]+=N
    return h

def main() -> None:
    if OUT.exists():
        raise SystemExit(f"refusing to overwrite existing archive {OUT}")
    OUT.mkdir(parents=True)
    p4=ENTRY/"prefix-full-rows2to4-v2"; p57=ENTRY/"prefix-full-rows5to7-v1"; root=ENTRY/"root-and-trees-v1"
    for d,expect in ((p4,1140),(p57,1140)):
        if len(list((d/"cores").glob("*.json"))) != expect: raise RuntimeError(f"incomplete {d}")
    s4=json.loads((p4/"summary.json").read_text()); s57=json.loads((p57/"summary.json").read_text()); sr=json.loads((root/"summary.json").read_text())
    if not (s4["status"]==s57["status"]==sr["status"]=="complete"): raise RuntimeError("incomplete source summary")

    # Copy exact, small source/summary artifacts.  The big per-core JSON trees
    # are represented below by canonical records plus source-byte ledgers.
    frozen=[]
    for src,rel in [
        (ENTRY/"s1_degree8_prefix.py","source/s1_degree8_prefix.py"),
        (ENTRY/"s1_degree8_stage5to7.py","source/s1_degree8_stage5to7.py"),
        (ENTRY/"s1_degree8_root_and_trees.py","source/s1_degree8_root_and_trees.py"),
        (ENTRY/"verify_s1_degree8_root_independent.py","source/verify_s1_degree8_root_independent.py"),
        (p4/"summary.json","frozen/production-stage4-summary.json"),
        (p57/"summary.json","frozen/production-stage57-summary.json"),
        (root/"summary.json","frozen/production-root-summary.json"),
        (ENTRY/"root-independent-v1.json","frozen/independent-root-summary.json"),
        (ALG/"s1_degree8_independent.py","dependencies/algebraic/s1_degree8_independent.py"),
        (ALG/"compare_full_records.py","dependencies/algebraic/compare_full_records.py"),
        (ALG/"full-v1/summary.json","dependencies/algebraic/full-v1-summary.json"),
        (ALG/"full-v1/comparison.json","dependencies/algebraic/full-v1-comparison.json"),
    ]: frozen.append(copy_frozen(src,rel))
    put_json("frozen/copied-artifacts.json", frozen)

    # Stage 4: canonical core list, compact count ledger, and all 60,134 rows.
    p4files=sorted((p4/"cores").glob("*.json"))
    p57files=sorted((p57/"cores").glob("*.json"))
    p4rows=[]; p4count=[]; cores=[]
    for f in p4files:
        x=json.loads(f.read_text()); core=x["binding"]["core"]; cores.append(core)
        rec=x["records"]
        p4rows.append({"core":core,"records":rec})
        p4count.append({"core":core,"records":len(rec),"source_name":f.name,"source_sha256":sha_path(f),"source_bytes":f.stat().st_size,"status":x["status"]})
    if len(set(tuple(x) for x in cores))!=1140 or sum(x["records"] for x in p4count)!=60134: raise RuntimeError("stage4 coverage error")
    put_json("data/core-list.json", {"cores":cores,"source_core_list_sha256":s4["core_list_sha256"],"count":len(cores)})
    put_gz_jsonl("data/stage4-core-counts.jsonl.gz",p4count)
    stage4parts=[]
    for start in range(0,len(p4rows),100):
        rel=f"data/stage4-records-{start//100:02d}.jsonl.gz"; p,raw=put_gz_jsonl(rel,p4rows[start:start+100]); stage4parts.append({"file":rel,"cores":len(p4rows[start:start+100]),"raw_sha256":raw,"records":sum(len(q["records"]) for q in p4rows[start:start+100])})
    put_json("data/stage4-parts.json",stage4parts)

    # Stage 5--7: every core count plus all 67 row-seven records.
    p57_by_core={}; p57count=[]; p57rows=[]; stage_mothers=[]
    for f in p57files:
        x=json.loads(f.read_text()); core=x["binding"]["core"]; p57_by_core[tuple(core)]=x
        rec=x["records"]; p57rows.append({"core":core,"records":rec})
        T,N,_ = core
        stage_mothers.extend([[T,N,*r[:6]] for r in rec])
        p57count.append({"core":core,"records":len(rec),"source_name":f.name,"source_sha256":sha_path(f),"source_bytes":f.stat().st_size,"status":x["status"]})
    if set(p57_by_core)!=set(tuple(x) for x in cores) or sum(x["records"] for x in p57count)!=67: raise RuntimeError("stage57 coverage error")
    put_gz_jsonl("data/stage57-core-counts.jsonl.gz",p57count)
    p,raw57=put_gz_jsonl("data/stage57-records.jsonl.gz",p57rows)
    put_json("data/stage57-records-index.json",{"file":"data/stage57-records.jsonl.gz","raw_sha256":raw57,"cores":1140,"records":67})

    # Exact production root output ledger, augmented (not substituted) with H.
    prod_by_row={}
    for f in sorted(root.glob("mother-*.json")):
        x=json.loads(f.read_text()); row=x["binding"]["row"]; key=tuple(row)
        if key in prod_by_row: raise RuntimeError("duplicate production mother")
        prod_by_row[key]={"source_name":f.name,"source_sha256":sha_path(f),"source_bytes":f.stat().st_size,"screen":x["screen"]}
    unique=[]; seen_mothers=set()
    for row in stage_mothers:
        if tuple(row) not in seen_mothers:
            seen_mothers.add(tuple(row)); unique.append(row)
    if len(unique)!=67 or set(map(tuple,unique)) != set(prod_by_row): raise RuntimeError("root-stage57 mapping error")
    mothers=[]
    for row in sorted(unique):
        pr=prod_by_row[tuple(row)]
        if pr["screen"]["status"]!="rejected": raise RuntimeError("unexpected root pass")
        mothers.append({"row":row,"H":h_from_row(row),"failure_code":pr["screen"]["code"],"production":pr})
    p,rawroot=put_gz_jsonl("data/root-mothers.jsonl.gz",mothers)
    put_json("data/root-mothers-index.json",{"file":"data/root-mothers.jsonl.gz","raw_sha256":rawroot,"mothers":67,"codes":dict(sorted(Counter(x["failure_code"] for x in mothers).items()))})

    # Byte-level source ledgers make external, un-copied raw per-core outputs
    # auditable without placing 2,280 redundant JSON files into a public tree.
    put_gz_jsonl("frozen/stage4-original-core-ledger.jsonl.gz",source_ledger(p4))
    put_gz_jsonl("frozen/stage57-original-core-ledger.jsonl.gz",source_ledger(p57))
    algcores=source_ledger(ALG/"full-v1")
    put_gz_jsonl("dependencies/algebraic/full-v1-core-ledger.jsonl.gz",algcores)

    status={
      "scope":"Compact archival index only; it validates stored bytes, coverage and cross-links, not the underlying mathematical enumeration.",
      "stage4":{"cores":1140,"records":60134,"summary_sha256":sha_path(p4/"summary.json")},
      "stage57":{"cores":1140,"records":67,"summary_sha256":sha_path(p57/"summary.json")},
      "root":{"mothers":67,"all_rejected":True,"summary_sha256":sha_path(root/"summary.json"),"independent_root_summary_sha256":sha_path(ENTRY/"root-independent-v1.json")},
      "algebraic_dependency":{"cores":1140,"summary_sha256":sha_path(ALG/"full-v1/summary.json"),"comparison_sha256":sha_path(ALG/"full-v1/comparison.json")},
    }
    put_json("archive-summary.json",status)

    readme="""# Round 055, degree-8 S={1}: compact archival index\n\nThis directory freezes the sources, summaries, canonical records, source-byte ledgers and independent-comparison bindings used by the necessary-prefix/root-screen work.  `verify_archive_index.py` checks archive hashes, gzip payloads, cardinalities, stage links, root mother/H/code links and the copied independent comparison summaries.  It is **not** a new complete mathematical enumeration or a substitute for the parent portable full verifier.\n\nThe 1,140 stage-4 and stage-5--7 original per-core JSON outputs are intentionally represented by canonical records plus source-byte ledgers: copying them would add redundant small-file noise.  All 60,134 row-4 records, all 67 row-7 records and all 67 root mothers (with exact `H=P*Pbar` and failure code) are retained.  Every gzip member is below 10 MB.\n\nHistorical disclosure: an early, unarchived independent-root draft used an incorrect `N*x^2*F` term. It disagreed immediately and was replaced before any passing certificate; it is not in this package.  The archived independent root verifier uses `F^2 + T*x*F + N*x^2`.  The historical stage-4 source has an integer-key JSON-resume binding issue (JSON makes those keys strings); it is frozen byte-for-byte here, and later portable verification must use a fresh run rather than trust that cache.\n\nRun only the archive/index check:\n\n```sh\npython3 verify_archive_index.py\n```\n\nThis package uses no symlinks and does not depend on the original absolute paths for verification.  Source locations and source hashes in the ledgers are provenance information.\n"""
    (OUT/"README.md").write_text(readme)
    # Include builders after all data; manifest excludes only itself by design.
    for src,rel in [(SELF/"build_package.py","source/build_package.py"),(SELF/"verify_archive_index.py","verify_archive_index.py")]: copy_frozen(src,rel)
    files=[]
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p.name!="package-manifest.json":
            if p.is_symlink(): raise RuntimeError(f"archive symlink {p}")
            files.append({"path":str(p.relative_to(OUT)),"sha256":sha_path(p),"bytes":p.stat().st_size})
    put_json("package-manifest.json",{"format":1,"files":files,"self_hash_excluded":True,"note":"The manifest is excluded from its own recursive hash list; every other regular file is listed."})
    print(json.dumps({"status":"complete","package":str(OUT),"files":len(files),"manifest_sha256":sha_path(OUT/"package-manifest.json")},sort_keys=True))

if __name__ == "__main__": main()
