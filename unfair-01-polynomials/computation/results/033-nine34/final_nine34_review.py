"""Read-only terminal review of theorem interfaces plus the final Fraction certificate."""
import hashlib,json,subprocess,sys
from pathlib import Path
SYNC=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/sync/unfair-01-polynomials')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 proof=SYNC/'proofs/033-nine-point-total34.md';root=SYNC/'computation/results/033-nine34';bind=json.loads((root/'parent-binding.json').read_text());all46=json.loads((root/'parent-all46-verification.json').read_text());cert=root/'positive/cached-recovered-candidate.json';d=json.loads(cert.read_text())
 text=proof.read_text();assert '33个次数对' in text and '严格在 \\((0,1)\\)' in text and '字节和覆盖绑定，不冒充再次运行' in text
 assert bind['prior_full_hash_receipt_sha256']=='ba8ee877c781ae9fc878c14a186cb83fa88a42a9bed1b93f24bbb02fdb025ccf' and len(bind['copied_files'])==68
 assert all46['status']=='all_46_residual_supports_exactly_certified' and all46['direct_count']==39 and all46['positive_count']==7 and len(all46['rows'])==46
 I,J=d['case']['I'],d['case']['J'];strict=set(d['strict_variables']);assert strict=={f'a{i}' for i in I[1:-1]}|{f'b{j}' for j in J[1:-1]};assert d['unit_variables']==['a0','a16','b0','b18']
 run=subprocess.run([sys.executable,str(SYNC/'computation/verify_positive_certificate.py'),str(cert)],text=True,capture_output=True,check=True);result=json.loads(run.stdout);assert result['status']=='verified'
 out={'status':'nine34_final_theorem_interface_review_passed','proof_sha256':sha(proof),'parent_binding_sha256':sha(root/'parent-binding.json'),'all46_sha256':sha(root/'parent-all46-verification.json'),'cached_certificate_sha256':sha(cert),'fraction_final_certificate':result,'notes':'Checked theorem wording, 33-pair interface, strict-variable boundary, and one final certificate. Hash bindings are recorded as bindings, not a mathematical re-run.'}
 p=Path(__file__).with_name('final-nine34-review.json');p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
