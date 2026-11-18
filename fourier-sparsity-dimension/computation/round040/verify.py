"""Direct Walsh checks of all recorded minimum bases and MM witnesses."""
import json
from pathlib import Path

here=Path(__file__).resolve().parent
def phase(i,j,n):
    return sum((((x>>i)&1)*((x>>j)&1))<<x for x in range(1<<n))
def check(base,generators,variables,expected_support):
    points=1<<variables
    masks=[sum(((x&i).bit_count()%2)<<x for x in range(points)) for i in range(points)]
    words=[0]
    support=0
    for g in generators:
        support|=g
        words += [g^w for w in words]
    assert support.bit_count()==expected_support and len(set(words))==len(words)
    assert all(abs(points-2*(base^w^m).bit_count())==1<<(variables//2) for w in words for m in masks)
census=json.loads((here/'census_4d_result.json').read_text())['rows']
assert len(census)==28
for row in census:
    for d,expected,count in [(4,12,1500),(5,14,5400),(6,15,720),(7,16,45)]:
        entry=row['dimensions'][str(d)]
        size,*generators=entry['minimum_and_basis']
        assert len(generators)==d and size==expected and not any(entry['hist'][:size]) and entry['hist'][size]==count
        check(row['seed'],generators,4,size)
mm=json.loads((here/'macbook_mm_e4_five_result.json').read_text())
base=phase(0,3,6)^phase(1,4,6)^phase(2,5,6)
for row in mm['extensions']:
    generators=[0x0101010101010101,0x1010101010101010,0x1111111100000000]
    generators.extend(int(w,16) for w in row['witness'])
    check(base,generators,6,row['joint_support'])
print('verified 112 four-variable bases and five six-variable extensions')
