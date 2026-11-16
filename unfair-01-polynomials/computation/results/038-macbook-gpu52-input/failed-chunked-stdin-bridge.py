import hashlib,json,os,subprocess,time
from pathlib import Path
assert os.path.ismount('/Volumes/data')
r=Path(__file__).resolve().parent
data=(r/'inputs.json').read_text()
wl='WriteString[$Output,"A30_STREAM_START\\n"];Flush[$Output];\ninputParts={};\n'
for offset in range(0,len(data),8192):
 wl+='AppendTo[inputParts,'+json.dumps(data[offset:offset+8192])+'];\n'
wl+='data=ImportString[StringJoin[inputParts],"RawJSON"];Clear[inputParts];\n'+(r/'search.wl').read_text()
cmd=['ssh','-S','/Volumes/data/CodexTemp/T2/unfair-01-polynomials/macbook-control.sock','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','macbook-math',"/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt"]
with (r/'stdout.jsonl').open('w') as out,(r/'stderr.txt').open('w') as err:
 p=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=out,stderr=err,text=True)
 p.stdin.write(wl);p.stdin.close();rc=p.wait()
(r/'exit.json').write_text(json.dumps({'returncode':rc,'input_sha256':hashlib.sha256(data.encode()).hexdigest(),'program_sha256':hashlib.sha256(wl.encode()).hexdigest(),'command':cmd,'finished_unix':time.time()},indent=2)+'\n')
raise SystemExit(rc)
