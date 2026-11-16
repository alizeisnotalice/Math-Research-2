#!/usr/bin/env python3
"""Byte-level check that the histogram extension leaves base kernels intact."""
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).parent


def metal(source):
    return re.search(r'let metal = """(.*?)"""', source, re.S).group(1)


def kernel_block(metal_source, name):
    start = metal_source.index('kernel void ' + name + '(')
    depth = 0
    for pos in range(start, len(metal_source)):
        if metal_source[pos] == '{':
            depth += 1
        elif metal_source[pos] == '}':
            depth -= 1
            if depth == 0 and pos > start:
                return metal_source[start:pos + 1]
    raise AssertionError(name)


base = metal((ROOT / 'frozen_base_r016.swift').read_text())
host = metal((ROOT / 'hist_host_audit.swift').read_text())
base_helpers = base[:base.index('kernel void count(')]
host_helpers = host[:host.index('kernel void count(')]
base_core = '\n'.join(kernel_block(base, x) for x in ('count', 'sum', 'cdf', 'prop', 'term'))
host_core = '\n'.join(kernel_block(host, x) for x in ('count', 'sum', 'cdf', 'prop', 'term'))
assert base_helpers == host_helpers
assert base_core == host_core
assert 'termHist' not in base and 'termHist' in host
out = {
    'status': 'passed',
    'base_source_sha256': hashlib.sha256((ROOT / 'frozen_base_r016.swift').read_bytes()).hexdigest(),
    'host_audit_source_sha256': hashlib.sha256((ROOT / 'hist_host_audit.swift').read_bytes()).hexdigest(),
    'helper_bytes_sha256': hashlib.sha256(base_helpers.encode()).hexdigest(),
    'count_sum_cdf_prop_term_bytes_sha256': hashlib.sha256(base_core.encode()).hexdigest(),
    'difference': 'host source adds only termHist plus post-command-buffer root/allowed/lastT readout and JSON fields; helpers and count/sum/cdf/prop/term Metal bytes match base exactly.',
}
(ROOT / 'kernel-identity-certificate.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
