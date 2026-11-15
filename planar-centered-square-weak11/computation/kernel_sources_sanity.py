import math, json, sys, hashlib, pathlib
root=pathlib.Path(__file__).resolve().parents[1]
files=["proofs/kernel-interface-audit.md","proofs/kernel-family-floor.md","references/kernel-interface-versions.json"]
out={"python":sys.version,"scope":"formula evaluation only; no C4 or LP rerun", "source_hashes":{f:hashlib.sha256((root/f).read_bytes()).hexdigest() for f in files}}
rows=[]
for alpha in (0.5,1.0,1.175,1.5,1.9):
 R=(2/alpha)**(1/(2-alpha))
 J=math.pi*(2-alpha)*R**2/(2*alpha*(R**(2-alpha)-1))
 rows.append({"alpha":alpha,"optimal_R":R,"isotropic_truncated_potential_halfmass":J,"same_expression":math.pi*R*R/2})
 assert abs(J-math.pi*R*R/2)<1e-12
out["isotropic_stable_formula_evaluations"]=rows
out["brownian_log_kernel_R_sqrt_e_halfmass"]=math.pi*math.e/2
print(json.dumps(out,indent=2))
