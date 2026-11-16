#!/usr/bin/env python3
"""Exact terminal-group overflow bounds for the r020 n<=16 contract."""
import json
import pathlib
from math import isqrt

n=16
group=256
# At r=N a live particle has one unoccupied vertex, so its allowed count a is 1.
coordinate=n-1
norm_max=2*coordinate*coordinate
f_per_particle=2*norm_max
g_abs_per_particle=norm_max
certificate={
    "n_max":n,"terminal_allowed_weight_max":1,"group_particles_max":group,
    "coordinate_abs_max":coordinate,"root_norm_sq_max":norm_max,
    "F_per_particle_max":f_per_particle,"G_abs_per_particle_max":g_abs_per_particle,
    "F_group_max":group*f_per_particle,"G_group_abs_max":group*g_abs_per_particle,
    "uint64_max":2**64-1,"int64_max":2**63-1,
    "F_uint64_safe":group*f_per_particle < 2**64,
    "G_int64_safe":group*g_abs_per_particle < 2**63,
}
assert certificate["F_uint64_safe"] and certificate["G_int64_safe"]
full=[]
for side in range(2,17):
    N=side*side
    # Exact consequence of the host guard P^2*N <= UInt64.max.
    pmax=isqrt((2**64-1)//N)
    fmax=pmax*4*(side-1)*(side-1)
    gmax=pmax*2*(side-1)*(side-1)
    full.append({"n":side,"N":N,"P_max":pmax,"F_total_max":fmax,"G_total_abs_max":gmax,
                 "F_uint64_safe":fmax<2**64,"G_int64_safe":gmax<2**63})
assert all(x["F_uint64_safe"] and x["G_int64_safe"] for x in full)
certificate["full_particle_terminal_bounds_under_P2N_guard"]=full
path=pathlib.Path(__file__).with_name("terminal-integer-bounds.json")
path.write_text(json.dumps(certificate,indent=2)+"\n")
print(json.dumps(certificate))
