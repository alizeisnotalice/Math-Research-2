import base64,gc,hashlib,json,os,subprocess,time,sys
all_data=json.loads(base64.b64decode(INPUT_B64));search=base64.b64decode(SEARCH_B64).decode();started=time.monotonic();completed=0
emit=lambda r:print(json.dumps(r,separators=(',',':')),flush=True)
emit(dict(event='controller_started',pid=os.getpid(),python=sys.version,tasks=len(all_data['tasks']),seconds_limit=14400,storage='RAM only; no remote task files'))
for task in all_data['tasks']:
 remaining=14400-(time.monotonic()-started)
 if remaining<=1:break
 data={**all_data,'tasks':[task]}
 wl='data=ImportString['+json.dumps(json.dumps(data,separators=(',',':')))+',"RawJSON"];\n'+search
 assert len(wl.encode())<245000
 p=subprocess.Popen(['/Applications/Wolfram.app/Contents/MacOS/WolframKernel','-noprompt','-run',wl],stdin=subprocess.DEVNULL)
 emit(dict(event='kernel_started',task=task['id'],pid=p.pid,program_sha256=hashlib.sha256(wl.encode()).hexdigest(),program_bytes=len(wl.encode())))
 timed_out=False
 try:rc=p.wait(timeout=remaining)
 except subprocess.TimeoutExpired:
  timed_out=True;p.terminate()
  try:rc=p.wait(timeout=5)
  except subprocess.TimeoutExpired:p.kill();rc=p.wait()
 emit(dict(event='kernel_reaped',task=task['id'],pid=p.pid,returncode=rc,timed_out=timed_out))
 del data,wl;gc.collect()
 if rc!=0 or timed_out:break
 completed+=1
emit(dict(event='controller_complete' if completed==len(all_data['tasks']) else 'controller_partial',completed_tasks=completed,tasks=len(all_data['tasks']),seconds=time.monotonic()-started,all_started_kernels_reaped=True))
