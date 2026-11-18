"""Check all disjoint work units and the closed fundamental-cut counts."""
from pathlib import Path

root=Path(__file__).resolve().parent
m=12
total=2*m*3**(m-1)
rows=[list(map(int,row.split())) for file in ('mini.txt','macbook.txt')
      for row in (root/file).read_text().splitlines()]
assert sorted(row[0] for row in rows)==list(range(108))
assert all(row[1]==len(range(row[0],total,108)) and row[2]+row[3]==row[1]*(2*m-1) for row in rows)
two=sum(row[2] for row in rows)
four=sum(row[3] for row in rows)
# Per tree, the broken segment always contributes one degree-one internal vertex.
# Each ordinary segment has such a vertex in two of its three tree choices.
assert 3*two==total*(3+2*(m-1))
assert two+four==total*(2*m-1)
print(dict(spanning_trees=total,fundamental_cuts=two+four,size_two=two,size_four=four))
