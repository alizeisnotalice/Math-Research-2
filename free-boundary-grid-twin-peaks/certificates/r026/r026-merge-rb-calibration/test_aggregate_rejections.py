#!/usr/bin/env python3
"""Negative input-audit checks: wrong B baseline and missing/reordered RB IDs."""
import json,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).parent
RAW=ROOT/'raw/n3-p4096-r64.jsonl'; SRC=ROOT/'source/tau2_merge_rb_calibration_r026.swift'
TOP=ROOT/'exact-top2-d4'
def rejected(raw,top):
 return subprocess.run([sys.executable,ROOT/'aggregate_calibration.py',raw,SRC,top],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).returncode!=0
with tempfile.TemporaryDirectory() as d:
 d=Path(d); lines=RAW.read_text().splitlines(); bad=d/'missing.jsonl'; kept=[x for x in lines if not (json.loads(x).get('kind')=='firstmerge_rb_calibration' and json.loads(x).get('replicate')==3)]; bad.write_text('\n'.join(kept)+'\n'); assert rejected(bad,TOP)
 wrong=d/'top'; wrong.mkdir(); (wrong/'summary-n3-n7.json').write_bytes((TOP/'summary-n3-n7.json').read_bytes()); x=json.loads((TOP/'n3.json').read_text());x['Btop2']='1';(wrong/'n3.json').write_text(json.dumps(x)); assert rejected(RAW,wrong)
print(json.dumps({'passed':True,'negative_cases':['missing RB record rejected','wrong exact B rejected']}))
