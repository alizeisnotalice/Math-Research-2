from pathlib import Path
import json,subprocess,hashlib,sys
r=Path(__file__).resolve().parent;tmp=Path('/Volumes/data/CodexTemp/T2/a022-r019-gpu-audit-reproduce');tmp.mkdir(parents=True,exist_ok=True)
subprocess.run(['swiftc','-O','-module-cache-path','/Volumes/data/CodexTemp/cache/swift-modules',str(r/'soa_frozen.swift'),'-o',str(tmp/'soa_frozen')],check=True)
old=json.loads((r/'root-frozen-62-audit.json').read_text());records=[]
sh=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
source_sha=sh(r/'soa_frozen.swift');verify_sha=sh(r/'verify_soa_frozen.py')
assert source_sha==old['source_sha256'] and verify_sha==old['verifier_sha256']
for z in old['records']:
 out=tmp/f"n{z['n']}-p{z['P']}-t{z['tau']}-l{z['L']}.jsonl"
 cmd=[str(tmp/'soa_frozen'),str(z['n']),str(z['P']),'1',str(z['tau']),str(z['tau']),str(out),'1','strat',str(z['L']),'trace']
 run=subprocess.run(cmd,capture_output=True,text=True);assert run.returncode==0,(cmd,run.stderr)
 ver=subprocess.run([sys.executable,str(r/'verify_soa_frozen.py'),str(out)+'.trace.json'],capture_output=True,text=True);assert ver.returncode==0,(cmd,ver.stderr)
 records.append({'n':z['n'],'P':z['P'],'tau':z['tau'],'L':z['L'],'command':cmd,'run_rc':0,'verify_rc':0,'output_sha256':sh(out),'trace_sha256':sh(Path(str(out)+'.trace.json')),'verification':json.loads(ver.stdout)})
result={'source_sha256':source_sha,'verifier_sha256':verify_sha,'binary_sha256':sh(tmp/'soa_frozen'),'records':records,'scope':'Root independently rebuilt frozen source, reran all 62 actual GPU cases, verified current exact trace hashes.'}
(tmp/'reproduced-62-audit.json').write_text(json.dumps(result,indent=2)+'\n');print('root actual GPU run and current verifier passed:',len(records))
