import pathlib,argparse,subprocess,sys,os,time,json
p=argparse.ArgumentParser();p.add_argument('--output',type=pathlib.Path,required=True);z=p.parse_args();assert str(z.output.resolve()).startswith('/Volumes/data/CodexTemp/');z.output.mkdir(exist_ok=False)
b=pathlib.Path(__file__).resolve().parent;records=[]
for name,src in [('block','block-encoding/audit_block_encoding_n3.py'),('coefficients','block-encoding/audit_coefficients.py'),('root-all-graphs','root-comparison/audit_all_graphs.py'),('boundary','boundary-bridge/audit_boundary_solechild.py'),('image','block-encoding/audit_image_n3.py')]:
 start=time.monotonic()
 with (z.output/(name+'.stdout')).open('w') as so,(z.output/(name+'.stderr')).open('w') as se:
  c=subprocess.Popen([sys.executable,'-B',str(b/src)],stdout=so,stderr=se);rc=c.wait()
 try:os.kill(c.pid,0);absent=False
 except ProcessLookupError:absent=True
 records.append(dict(name=name,pid=c.pid,exit_code=rc,waited_reaped=True,absent_after=absent,wall_seconds=time.monotonic()-start));(z.output/'receipt.json').write_text(json.dumps(dict(python=sys.version,records=records),indent=2)+'\n');assert rc==0 and absent,name
print('PASS',len(records),'checks; no production started')
