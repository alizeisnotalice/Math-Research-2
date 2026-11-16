#!/usr/bin/env python3
"""Integer checks for the r019 canonical root-image ledgers."""
import json
from fractions import Fraction
from pathlib import Path

expected={
 3:(29568,143808,97618),
 4:(35704394880,779299174080,389050636544),
 5:(88366736882654697600,7063261481362958280000,2512023954782991730592),
}
rows=[]
for n,(a1,a2,image) in expected.items():
 d=json.loads(Path(f'root-image-n{n}-r019.json').read_text())
 assert (int(d['A1']),int(d['A2_from_full_first_merge_scan']),int(d['canonical_root_image_sum']))==(a1,a2,image)
 sizes=d['by_canonical_B_size']; ranks=d['by_first_merge_rank']
 assert sum(int(x['W']) for x in sizes)==a2
 assert sum(int(x['image']) for x in sizes)==image
 assert sum(int(x['W']) for x in ranks)==a2
 assert sum(int(x['image']) for x in ranks)==image
 rows.append({'n':n,'E1J':str(Fraction(image,a1)),'E1J_over_N':str(Fraction(image,n*n*a1)),'E2_qB':str(Fraction(image,a2))})
out={'status':'passed','rows':rows}
lobes=json.loads(Path('scattered-lobes-n3.json').read_text())
assert lobes['status']=='passed' and lobes['A1']==29568 and lobes['events']==0
assert lobes['weighted_source_sum_numerator']==0 and lobes['max_weighted_output_load_numerator']==0
out['scattered_lobes_n3']={'events':0,'status':'passed'}
subsets=json.loads(Path('scattered-subsets-n3.json').read_text())
assert subsets['status']=='passed' and subsets['A1']==29568 and subsets['invalid_outputs']==0
assert (subsets['weighted_sum_num'],subsets['weighted_sum_den'])==(534178,7)
assert Fraction(subsets['max_load_num'],subsets['max_load_den'])<=1
out['scattered_subsets_n3']={'events':subsets['events'],'weighted_sum':'534178/7','status':'passed'}
Path('verify-r019.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
