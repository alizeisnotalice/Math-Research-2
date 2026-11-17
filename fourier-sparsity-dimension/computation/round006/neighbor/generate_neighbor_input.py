#!/usr/bin/env python3
"""Generate the audited corrected-results row-1 constants for neighbour search."""
import argparse, hashlib, json
from pathlib import Path

p=argparse.ArgumentParser(); p.add_argument('--audit',type=Path,required=True); p.add_argument('--partitions',type=Path,required=True); p.add_argument('--out',type=Path,required=True); a=p.parse_args()
ra,rp=a.audit.read_bytes(),a.partitions.read_bytes(); audit,parts=json.loads(ra),json.loads(rp)
row=next(r for r in audit['rows'] if r['source_file'].endswith('corrected-results.json') and r['row_index']==1)
part=next(r for r in parts['rows'] if r['source_file']==row['source_file'] and r['row_index']==1)
if len(row['relative_generators'])!=8 or len(part['partition_masks_hex'])!=8: raise AssertionError('shape')
fmt=lambda xs: ','.join(f'0x{int(x):016x}ULL' for x in xs)
blocks=[int(x,16) for x in part['partition_masks_hex']]
a.out.write_text(f'''// generated; audit_sha256={hashlib.sha256(ra).hexdigest()} partitions_sha256={hashlib.sha256(rp).hexdigest()}\n#pragma once\n#include <array>\n#include <cstdint>\nconstexpr std::array<std::uint64_t,8> INPUT_REL={{ {fmt(row['relative_generators'])} }};\nconstexpr std::array<std::uint64_t,8> INPUT_H={{ {fmt(blocks)} }};\nconstexpr std::array<int,8> INPUT_PERM={{ {','.join(map(str,row['permutation']))} }};\nconstexpr int INPUT_H_BITS={row['h_bits']};\n''')
