#!/usr/bin/env python3
"""Exact cross-coordinate S64 audit and verified all-q independent search."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
from time import monotonic

P=163
ROOT=Path(__file__).resolve().parent
REFERENCE=ROOT/"reference_source.cpp"
FIXED_ONES=(31,32,33,34,35,38,58)
FIXED_ZEROS=(60,61,62,63)


def encoded(x):
    return json.dumps(x,sort_keys=True,separators=(",",":")).encode()+b"\n"


def digest(data):
    return sha256(data).hexdigest()


def file_digest(p):
    return digest(p.read_bytes())


def trim(A):
    while A and A[-1]==0:
        A.pop()
    return A


def divide(A,B,p=P):
    work=trim([x%p for x in A] if p is not None else list(map(Fraction,A)))
    quotient=[0]*max(0,len(work)-len(B)+1)
    inverse=pow(B[-1],-1,p) if p is not None else 1/Fraction(B[-1])
    while len(work)>=len(B):
        k=len(work)-len(B)
        c=work[-1]*inverse
        if p is not None:c%=p
        quotient[k]=c
        for j,b in enumerate(B):
            value=work[k+j]-c*b
            work[k+j]=value%p if p is not None else value
        trim(work)
    return trim(quotient),work


def squarefree(A,p=P):
    B=[i*A[i] for i in range(1,len(A))]
    if p is not None:
        A=[x%p for x in A]
        B=[x%p for x in B]
    A,B=trim(A),trim(B)
    while B:
        _,r=divide(A,B,p)
        A,B=B,r
    return len(A)==1


def J(q):
    return [q,117,115,68,106,121,6,1]


def multiply(A,B,p=P):
    out=[0]*(len(A)+len(B)-1)
    for i,a in enumerate(A):
        for j,b in enumerate(B):
            out[i+j]=(out[i+j]+a*b)%p
    return out


def old_columns(degree,q):
    j2=multiply(J(q),J(q))
    rows=[]
    for i in range(degree+1):
        _,r=divide([0]*(degree-i)+[(-1)**i],j2)
        rows.append(r+[0]*(14-len(r)))
    return rows


def transformed(combo,columns):
    out=[0]*14
    for j,c in enumerate(combo):
        for k,x in enumerate(columns[j]):
            out[k]+=c*x
    return [x%P for x in out]


def wanted_indices(degree,mode,fixture=None):
    fixed={} if fixture is not None else {**{j:int(j!=7) for j in range(12)},degree:1}
    if mode=="tight":
        fixed.update(dict.fromkeys(FIXED_ONES,1))
        fixed.update(dict.fromkeys(FIXED_ZEROS,0))
    return fixed,[i for i in range(12,degree) if i not in fixed]


def verify_result(result,degree,q,mode,fixture=None):
    assert result["degree"]==degree and result["q"]==q and result["mode"]==mode
    fixed,free=wanted_indices(degree,mode,fixture)
    assert result["positions"]==free
    w=len(free)//2
    assert result["stored_patterns"]==1<<w and result["scanned_patterns"]==1<<(len(free)-w)
    expected_mask=sum(1<<(i-12) for i in fixed if 12<=i<degree and fixed[i]==1)
    assert result["fixed_mask"]==expected_mask
    cols=old_columns(degree,q)
    ordinary=result["column_remainders"]
    modulus=result["modulus"]
    assert len(ordinary)==degree+1
    assert all(len(v)==14 and all(0<=r<P for r in v) for v in ordinary)
    assert len(modulus)==15 and modulus[-1]==1 and all(0<=r<P for r in modulus)
    if degree>=14:
        assert transformed(modulus,cols)==[0]*14
    else:
        extended=old_columns(14,q)
        assert transformed(modulus,extended)==[0]*14
    for i,row in enumerate(ordinary):
        assert transformed(row,cols)==cols[i], (q,i)
    if fixture is None:
        target=[-fixed.get(i,0) for i in range(degree+1)]
    else:
        target=[(fixture>>(i-12))&1 if 12<=i<degree else 0 for i in range(degree+1)]
    assert transformed(result["target_remainder"],cols)==transformed(target,cols)
    projected=result["projected_masks"]
    assert len(projected)==result["projected_pairs"] and sorted(set(projected))==projected
    matches=result["matches"]
    independently_matched=[]
    survivors=[]
    axes=(0,3,6,9,12)
    for mask in projected:
        assert 0<=mask<1<<(degree-12)
        coeff=[(mask>>(i-12))&1 if 12<=i<degree else fixed.get(i,0) for i in range(degree+1)]
        assert all(coeff[i]==z for i,z in fixed.items())
        diff=[a-b for a,b in zip(coeff,target)] if fixture is not None else coeff
        assert all(sum(x*ordinary[i][j] for i,x in enumerate(diff))%P==0 for j in axes)
        if transformed(diff,cols)!=[0]*14:
            continue
        independently_matched.append(mask)
        if mode=="tight":
            quotient,R=divide(coeff,modulus)
            assert not R
            c_sqfree=any(squarefree(coeff,p) for p in (2,3,5,7,11,13,17,19)) or squarefree(coeff,None)
            if squarefree(quotient) and c_sqfree:survivors.append(mask)
    assert matches==independently_matched
    return dict(projected_mask_sha256=digest(encoded(projected)),column_sha256=digest(encoded(ordinary)),inverse_column_sha256=digest(encoded(cols)),
                result_sha256=digest(encoded(result)),max_square_degree7_char0_squarefree=survivors)


def small_case(bin_path,reference,degree,q,fixture):
    newargs=[str(bin_path),str(degree),str(q)]
    oldargs=[str(reference),str(degree),str(q)]
    if fixture is not None:
        newargs.append(f"--fixture={fixture}")
        oldargs.append(str(fixture))
    result=json.loads(subprocess.check_output(newargs))
    old=json.loads(subprocess.check_output(oldargs))
    summary=verify_result(result,degree,q,"fixture" if fixture is not None else "same_range",fixture)
    j2=multiply(J(q),J(q))
    fixed,free=wanted_indices(degree,"same_range",fixture)
    expected=[]
    for mask in range(1<<len(free)):
        bits=[((mask>>(i-12))&1) if i in free else fixed.get(i,0) for i in range(degree+1)]
        if fixture is not None:
            for i in free:bits[i]-=(fixture>>(i-12))&1
        _,r=divide([((-1)**i)*bits[i] for i in range(degree,-1,-1)],j2)
        if not r:expected.append(mask)
    assert result["matches"]==sorted(old["matches"])==expected
    if fixture is not None:assert fixture in expected
    return dict(degree=degree,q=q,fixture=fixture,matches=expected,bits=len(free),**summary)


def run_one(binary,q,path):
    begin=monotonic()
    command=["/usr/bin/time","-l",str(binary),"64",str(q),"--tight"]
    with subprocess.Popen(command,cwd=path.parent,stdout=subprocess.PIPE,stderr=subprocess.PIPE) as child:
        stdout,stderr=child.communicate(timeout=180)
        pid,code=child.pid,child.returncode
    assert code==0, (q,code,stderr)
    result=json.loads(stdout)
    checked=verify_result(result,64,q,"tight")
    rss=re.search(rb"(\d+)\s+maximum resident set size",stderr)
    assert rss, stderr
    output=dict(q=q,command=command,transport_pid=pid,exit_code=code,
                elapsed_seconds=round(monotonic()-begin,6),peak_rss_bytes=int(rss[1]),
                time_stderr=stderr.decode(),report=result,verification=checked)
    temporary=path.with_suffix(".partial")
    with temporary.open("xb") as f:
        f.write(encoded(output))
        f.flush()
        os.fsync(f.fileno())
    temporary.replace(path)
    return output


def main():
    if sys.flags.optimize:
        raise SystemExit("Python assertions must be enabled.")
    ap=argparse.ArgumentParser()
    ap.add_argument("--binary",type=Path,default=ROOT/"ordinary_remainder")
    ap.add_argument("--reference",type=Path,default=ROOT/"reference_small")
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    ap.add_argument("--smoke",action="store_true")
    a=ap.parse_args()
    out=a.out.resolve()
    if a.workers<1:ap.error("workers must be positive")
    out.mkdir(exist_ok=False)
    binary,reference=a.binary.resolve(),a.reference.resolve()
    q_values=list(range(1,P))
    assert all(squarefree(J(q)) for q in q_values)
    smoke=[small_case(binary,reference,*sample) for sample in
           ((13,1,None),(15,37,3),(18,162,None),(20,76,165),(23,110,1051),(24,2,None))]
    (out/"smoke.json").write_bytes(encoded(dict(status="verified",cases=smoke,brute_checks=sum(1<<s["bits"] for s in smoke),
                                                  reference_source_sha256=file_digest(REFERENCE))))
    if a.smoke:
        print(json.dumps(dict(status="verified",cases=len(smoke))))
        return
    started=monotonic()
    manifest=dict(degree=64,p=P,q=q_values,free_bits=41,reference_source_sha256=file_digest(REFERENCE),
                  source_sha256=file_digest(ROOT/"ordinary_remainder.cpp"),runner_sha256=file_digest(Path(__file__)),
                  binary_sha256=file_digest(binary),platform=platform.platform(),workers=a.workers)
    (out/"manifest.json").write_bytes(encoded(manifest))
    checks=out/"by-q"
    checks.mkdir()
    output=[]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures={q:pool.submit(run_one,binary,q,checks/f"q{q:03d}.json") for q in q_values}
        for q,future in futures.items():
            x=future.result()
            output.append(x)
            if q%10==0 or q==162:
                print(json.dumps(dict(completed_ordered_q=q,saved=len(list(checks.glob("*.json"))))),
                      flush=True)
    pids=[x["transport_pid"] for x in output]
    for pid in pids:
        try:
            os.kill(pid,0)
        except ProcessLookupError:
            pass
        else:
            raise RuntimeError(f"PID {pid} still present")
    assert [x["q"] for x in output]==q_values
    entries=[]
    for x in output:
        q=x["q"]
        r=x["report"]
        entries.append(dict(q=q,stored_patterns=r["stored_patterns"],scanned_patterns=r["scanned_patterns"],
                            projected_pairs=r["projected_pairs"],matches=r["matches"],**x["verification"],
                            checkpoint_sha256=file_digest(checks/f"q{q:03d}.json")))
    (out/"all-q.jsonl").write_bytes(b"".join(map(encoded,entries)))
    summary=dict(status="complete",parameters=len(output),free_bits=41,workers=a.workers,
                 total_stored=sum(r["stored_patterns"] for r in entries),
                 total_scanned=sum(r["scanned_patterns"] for r in entries),
                 exact_projection_pairs=sum(r["projected_pairs"] for r in entries),
                 full_matches=sum(len(r["matches"]) for r in entries),
                 full_scope_matches=sum(len(r["max_square_degree7_char0_squarefree"]) for r in entries),
                 total_cross_coordinate_columns=65*len(entries),
                 max_rss_per_child=max(x["peak_rss_bytes"] for x in output),
                 wall_seconds=round(monotonic()-started,6),tracked_processes_exited=len(pids),
                 sorted_q_index_sha256=file_digest(out/"all-q.jsonl"),
                 manifest_sha256=file_digest(out/"manifest.json"))
    (out/"complete.json").write_bytes(encoded(summary))
    print(json.dumps(summary,sort_keys=True,indent=2))


if __name__=="__main__":
    main()
