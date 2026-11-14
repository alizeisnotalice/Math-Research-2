#!/usr/bin/env python3
"""Transfer one millionth of row 0's mass from row 1 to row 0 of R12."""
import json,sys
from fractions import Fraction as F
from pathlib import Path

data=json.loads(Path(sys.argv[1]).read_text())
amount=F(data['pts'][0][2])/1000000
data['pts'][0][2]=str(F(data['pts'][0][2])+amount)
data['pts'][1][2]=str(F(data['pts'][1][2])-amount)
Path(sys.argv[2]).write_text(json.dumps(data,indent=2)+'\n')
