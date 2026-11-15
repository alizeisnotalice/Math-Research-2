#!/usr/bin/env python3
"""Restore archived search records to a new directory without overwriting evidence."""
import argparse, gzip, hashlib, json, os, shutil
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--archive',required=True,type=Path,help='topic root containing computation/')
p.add_argument('--output',required=True,type=Path,help='new, nonexisting output root')
a=p.parse_args(); archive=a.archive.resolve(); root=a.output.resolve()
if not os.path.ismount('/Volumes/data'): p.error('/Volumes/data is not mounted')
if not root.is_relative_to(Path('/Volumes/data/CodexTemp/T2')): p.error('output must be inside /Volumes/data/CodexTemp/T2')
if root.exists(): p.error('output root already exists; choose a new path')
comp=archive/'computation'; root.mkdir(parents=True)
records=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def copy(src,dst,origin):
 dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);records.append({'file':str(dst.relative_to(root)),'source':str(src.relative_to(archive)),'status':origin,'sha256':sha(dst)})
def unpack(src,dst):
 dst.parent.mkdir(parents=True,exist_ok=True)
 with gzip.open(src,'rb') as f,dst.open('xb') as out:shutil.copyfileobj(f,out)
 records.append({'file':str(dst.relative_to(root)),'source':str(src.relative_to(archive)),'status':'decompressed original bytes','sha256':sha(dst)})
copy(comp/'inputs/r008_final_0611.json',root/'round008/joint-groups/batch2/candidate_0611.json','archived final input bytes')
for number,sub,ident in [(25,'fixedmass-position','0110'),(26,'quadrant-groups','0360'),(27,'quadrant-joint','0079'),(28,'quadrant-joint','0026'),(29,'quadrant-joint','0298')]:
 folder=comp/'results'/f'round{number:03d}';dest=root/f'round{number:03d}'/sub
 for name in ('manifest.json','results.json'):unpack(folder/(name+'.gz'),dest/name)
 for name in ('summary.json','completed.jsonl'):
  if (folder/name).is_file():copy(folder/name,dest/name,'archived original bytes')
 name=f'r{number:03d}_final_{ident}.json';copy(comp/'inputs'/name,dest/f'candidate_{ident}.json','archived final input bytes')
 if not (dest/'completed.jsonl').exists():
  rows=json.loads((dest/'results.json').read_text());out=dest/'completed.reconstructed.jsonl'
  out.write_text(''.join(json.dumps(row,separators=(',',':'))+'\n' for row in rows))
  records.append({'file':str(out.relative_to(root)),'source':str((dest/'results.json').relative_to(root)),'status':'reconstructed sorted results; not original completion order or original bytes','sha256':sha(out)})
(root/'RESTORE-MANIFEST.json').write_text(json.dumps({'schema':'pcw-rounds-025-029-restoration-v1','files':records},indent=2)+'\n')
print(json.dumps({'output':str(root),'files':len(records),'manifest':str(root/'RESTORE-MANIFEST.json')},indent=2))
