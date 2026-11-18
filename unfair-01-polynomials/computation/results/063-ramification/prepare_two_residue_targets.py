"""Choose degree-56 Metal targets from words with two verified repeated residues."""

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import sympy as sp

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "062-residue-sieve"))

from prepare_gap7_degree56 import X, H, PRIMES, DEGREE, FREQUENCIES, RADIUS, factor_period, root_packet


HERE = Path(__file__).resolve().parent
MASKS = HERE / "ramified-masks-two-residues.txt"
SEED = 20261118057


def decode(mask):
    assert 0 <= mask < 1 << 44
    return [1]*7+[0]+[1]*4+[(mask >> bit) & 1 for bit in range(44)]+[1]


def value(c, root, modulus, derivative=False):
    if derivative:
        return sum(k*bit*pow(root,k-1,modulus) for k,bit in enumerate(c) if k) % modulus
    return sum(bit*pow(root,k,modulus) for k,bit in enumerate(c)) % modulus


def main():
    assert __debug__ and Path("/Volumes/data").is_mount()
    raw = MASKS.read_bytes()
    masks = list(map(int,raw.splitlines()))
    assert masks and len(masks) == len(set(masks))
    for mask in masks:
        c = decode(mask)
        assert value(c,-2,163) == value(c,-2,163,True) == 0
        assert value(c,2,1279) == value(c,2,1279,True) == 0
    rng = np.random.default_rng(SEED)
    rng.shuffle(masks)
    periods = {p: factor_period(H,p)[0] for p in PRIMES}
    points = RADIUS*np.exp(2j*np.pi*np.arange(FREQUENCIES)/FREQUENCIES)
    tasks, rejected = [], []
    for rank, mask in enumerate(masks):
        c = decode(mask)
        polynomial = sp.Poly.from_list(c[::-1], X, domain=sp.ZZ)
        if not polynomial.is_sqf:
            rejected.append({"mask":mask,"kind":"repeated_in_Q"})
            continue
        local = {}
        for p, h_period in periods.items():
            period, degrees = factor_period(c[::-1], p)
            if period % h_period:
                rejected.append({"mask":mask,"kind":"period","prime":p,"degrees":degrees})
                break
            local[p] = {"period":period,"degrees":degrees}
        else:
            result, reason = root_packet(np.array(c),points)
            if reason:
                rejected.append({"mask":mask,"kind":reason})
                continue
            for p in (163,1279):
                fp=sp.Poly(polynomial,modulus=p)
                assert sp.gcd(fp,fp.diff()).degree()>0
            tasks.append({"id":len(tasks),"mask":mask,"periods":local,**result})
            if len(tasks)==12:
                break
    print(json.dumps({"selected":len(tasks),"rejected":dict(Counter(r["kind"] for r in rejected))}),flush=True)
    assert tasks
    k=np.arange(DEGREE+1)
    phase=2*np.pi*np.outer(np.arange(FREQUENCIES),k)/FREQUENCIES
    packet={"degree":DEGREE,"L":FREQUENCIES,"seed":SEED,"batch":131072,"tasks":tasks,
            "DC":(np.cos(phase)/FREQUENCIES/(RADIUS**k)).tolist(),
            "DS":(np.sin(phase)/FREQUENCIES/(RADIUS**k)).tolist(),
            "h_reduction_periods":periods,"source_mask_list_sha256":sha256(raw).hexdigest(),
            "residue_roots":{"163":-2,"1279":2},"matched_words":len(masks),
            "attempted_draws":rank+1,"rejected_draws":rejected,
            "scope":"Selected numerical targets with the C7 gap and repeated roots at both ramified primes. "
                    "The rational squarefree condition is checked before numerical root allocation."}
    filename=HERE/"inputs-two-residue-degree56.json"
    content=json.dumps(packet,sort_keys=True,separators=(",",":"))+"\n"
    with filename.open("x") as f:
        f.write(content)
    result={"status":"verified","packet_sha256":sha256(content.encode()).hexdigest(),"candidates":len(tasks),
            "masks_total":sum(t["masks"] for t in tasks),"screened":rank+1,
            "generator_masks_verified":len(masks),"rejections":dict(Counter(r["kind"] for r in rejected))}
    with (HERE/"sieve-two-residue-check.json").open("x") as f:
        json.dump(result,f,indent=2)
        f.write("\n")
    print(json.dumps(result))


if __name__=="__main__":
    main()
