#!/Users/zzh/.venvs/main/bin/python
"""Independent exact free-prefix search using C=H*(GQ), not A*B."""
from __future__ import annotations
import gzip, hashlib, json, platform, sys, time
from math import isqrt
from pathlib import Path

SOURCE_SHA = "47961d58e3bbed6bff168fa19b3ec4dd5028945f5948e0ffebf79634de56d3a7"
MAX_DEPTH, STATE_CAP = 64, 100_000

def add(x, y): return (x[0]+y[0], x[1]+y[1])

def floor_alpha(z, T, delta):
    p, q = z
    if q == 0: return p
    u, r = 2*p+q*T, isqrt(q*q*delta)
    return (u+r)//2 if q > 0 else (u-r-1)//2

def options(seq, P, k, T, delta):
    # Includes P_k*seq[0] when k<=5.
    base = (0, 0)
    for j in range(1, min(k, 5)+1):
        base = add(base, (P[j][0]*seq[k-j], P[j][1]*seq[k-j]))
    f = floor_alpha(base, T, delta)
    return (-f, 1-f) if base[1] == 0 else (-f,)

def hcoeff(T, N, c, d):
    return (1,0,T,2*c,2*d+N,2+T*c,c*c+T*d,2*c*d+T,d*d+2*c,2*d,1)

def c_coefficient(G, Q, h, k):
    return sum(h[j]*sum(G[i]*Q[k-j-i] for i in range(k-j+1))
               for j in range(min(k, 10)+1))

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""): h.update(block)
    return h.hexdigest()

def save_gz(path, value):
    with gzip.open(path, "wt", encoding="utf-8") as f:
        json.dump(value, f, separators=(",", ":")); f.write("\n")

def run_one(index, core, output, checkpoint):
    T,N,c,d = core
    delta = T*T-4*N
    assert delta > 0 and isqrt(delta)**2 != delta
    P = ((1,0),(0,0),(0,1),(c,0),(d,0),(1,0))
    Pb = ((1,0),(0,0),(T,-1),(c,0),(d,0),(1,0))
    h, root = hcoeff(T,N,c,d), output/f"core-{index:03d}"
    root.mkdir(exist_ok=True)
    states, layers, status, death = [((1,),(1,))], [], "unknown_at_depth", None
    for k in range(1, MAX_DEPTH+1):
        accepted, rejected, tried, capped = [], [], 0, False
        for parent,(G,Q) in enumerate(states):
            for gk in options(G,P,k,T,delta):
                for qk in options(Q,Pb,k,T,delta):
                    tried += 1
                    G2,Q2 = G+(gk,),Q+(qk,)
                    ck = c_coefficient(G2,Q2,h,k)
                    if ck in (0,1):
                        accepted.append((G2,Q2))
                        if len(accepted) > STATE_CAP: capped=True; break
                    else: rejected.append([parent,gk,qk,ck])
                if capped: break
            if capped: break
        layer_file=root/f"layer-{k:02d}.json.gz"
        save_gz(layer_file,{"core":core,"depth":k,"tried":tried,
            "accepted":[{"G":list(G),"Q":list(Q)} for G,Q in accepted]})
        layers.append({"depth":k,"tried":tried,"accepted":len(accepted),
            "states_sha256":sha(layer_file),"capped":capped})
        checkpoint[str(index)]={"core":core,"layers":layers,"status":status}
        (output/"checkpoint.json").write_text(json.dumps(checkpoint,separators=(",",":"))+"\n")
        if capped: status="unknown_state_cap"; break
        if not accepted:
            death_file=root/f"death-{k:02d}.json.gz"
            save_gz(death_file,{"core":core,"depth":k,
                "previous_states":[{"G":list(G),"Q":list(Q)} for G,Q in states],
                "rejected":rejected})
            death={"depth":k,"tried":tried,"witness_sha256":sha(death_file)}
            status="dead_exhaustive_prefix"; break
        states=accepted
    result={"index":index,"core":core,"delta":delta,"status":status,"layers":layers,"death":death}
    (root/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    checkpoint[str(index)]=result
    (output/"checkpoint.json").write_text(json.dumps(checkpoint,separators=(",",":"))+"\n")
    return result

def main():
    here=Path(__file__).resolve().parent
    source=here.parent/"core5_survivors.json"
    assert sha(source)==SOURCE_SHA
    cores=json.loads(source.read_text())
    output=here/"hgq-depth64-cap100000"; output.mkdir(exist_ok=True)
    checkpoint={}; start=time.monotonic()
    results=[run_one(i,core,output,checkpoint) for i,core in enumerate(cores)]
    counts={}
    for row in results: counts[row["status"]]=counts.get(row["status"],0)+1
    summary={"status":"completed_independent_HGQ_free_prefix","source_sha256":SOURCE_SHA,
        "cores":len(cores),"max_depth":MAX_DEPTH,"state_cap":STATE_CAP,
        "status_counts":counts,"results":results,"seconds":time.monotonic()-start,
        "runtime":{"python":sys.version,"platform":platform.platform()},
        "scope":"C=H(GQ) includes h0*(Gk*Q0+G0*Qk). No tail is imposed. Dead prefixes exclude arbitrary finite degrees for that core; living/capped prefixes are unknown."}
    target=output/"summary.json"; target.write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps({"status_counts":counts,"seconds":summary["seconds"],"summary_sha256":sha(target)},indent=2))
if __name__=="__main__": main()
