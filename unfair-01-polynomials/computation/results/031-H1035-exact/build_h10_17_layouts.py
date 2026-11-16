"""Exact canonical list after the audited two-fiber trace reduction."""
import hashlib,json
from pathlib import Path
def main():
 # bit i=1 means the i-th conjugate mother pair takes its 5-fiber; bit 0=0
 # fixes the A/B exchange representative in the e=1 stratum.
 layouts=[{'e':0,'p5_mask':0,'kind':'constant'}]
 layouts += [{'e':1,'p5_mask':mask,'kind':'full_fibers'} for mask in range(32) if not(mask&1)]
 assert len(layouts)==17 and len({(x['e'],x['p5_mask']) for x in layouts})==17
 out={'status':'H10_two_fiber_17_canonical_layouts_built','layouts':layouts,'count':len(layouts),'scope':'Necessary layouts conditional only on the audited trace and root-product lemmas; no coefficient test.'}
 p=Path(__file__).with_name('H10-17-layouts.json');p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'count':17,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}))
if __name__=='__main__':main()
