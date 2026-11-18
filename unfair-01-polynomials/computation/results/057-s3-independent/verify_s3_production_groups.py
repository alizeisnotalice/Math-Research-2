#!/usr/bin/env python3
"""Compare production groups against frozen per-core tail shards, as sets."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def normalize_bounds(L, dlo, dhi, elo, ehi, zlo, zhi):
    slope = abs(L)
    if slope == 0:
        return dlo, dhi, max(elo,zlo), min(ehi,zhi)
    original_negative = L < 0
    if original_negative:
        dlo, dhi = -dhi,-dlo
    dlo = max(dlo, -((zhi-elo)//slope))
    dhi = min(dhi, (ehi-zlo)//slope)
    assert dlo <= dhi
    elo, ehi = max(elo,zlo+slope*dlo), min(ehi,zhi+slope*dhi)
    assert elo <= ehi
    if original_negative:
        dlo, dhi = -dhi,-dlo
    return dlo, dhi, elo, ehi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("production_groups", type=Path)
    ap.add_argument("tail_manifest", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    assert not a.out.exists()
    manifest = a.tail_manifest.read_bytes()
    metadata = [json.loads(line) for line in manifest.splitlines()]
    source_hash = hashlib.sha256()
    paths = 0
    source_bytes = 0
    groups = mothers = 0
    abs_z_extreme = 0
    with a.production_groups.open("rb") as production:
        for entry in metadata:
            if entry["path"] is None:
                assert entry["groups"] == entry["mothers"] == 0
                continue
            paths += 1
            path = a.tail_manifest.parent/entry["path"]
            payload = path.read_bytes()
            assert len(payload) == entry["bytes"] and hashlib.sha256(payload).hexdigest() == entry["sha256"]
            text = gzip.decompress(payload)
            assert len(text) == entry["raw_bytes"] and hashlib.sha256(text).hexdigest() == entry["raw_sha256"]
            group_count = mother_count = 0
            for line in text.splitlines():
                L,M,T,N,D,b,c,dlo,dhi,elo,ehi,zlo,zhi,count = json.loads(line)
                assert [L,M,T,N,D] == entry["core"]
                raw = production.readline()
                source_hash.update(raw)
                source_bytes += len(raw)
                values = tuple(map(int,raw.split()))
                assert values[:6] == (L,T,N,D,b,c) and values[-3:] == (zlo,zhi,count)
                assert normalize_bounds(L,*values[6:-1]) == (dlo,dhi,elo,ehi), (L,T,N,D,b,c)
                group_count += 1
                mother_count += count
                abs_z_extreme = max(abs_z_extreme,abs(zlo),abs(zhi))
            assert (group_count,mother_count) == (entry["groups"],entry["mothers"])
            groups += group_count
            mothers += mother_count
        assert production.readline() == b""
    out = dict(status="verified", cores=len(metadata), shards=paths, groups=groups, mothers=mothers,
               abs_z_extreme=abs_z_extreme, production_groups_bytes=source_bytes,
               production_groups_sha256=source_hash.hexdigest(),
               tail_manifest_sha256=hashlib.sha256(manifest).hexdigest())
    a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps(out,sort_keys=True))


if __name__ == "__main__":
    main()
