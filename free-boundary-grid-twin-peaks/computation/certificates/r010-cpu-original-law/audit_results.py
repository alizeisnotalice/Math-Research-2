#!/usr/bin/env python3
"""Compare CPU forced-rank SIS aggregates with exact small-grid ledgers."""
import json
import sys
from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 50
root = Path(__file__).parent
exact = {
    2: (Decimal(16), Decimal(8), Decimal(1)/2),
    3: (Decimal(29568), Decimal(143808), Decimal(1100)/Decimal(2247)),
    4: (Decimal(35704394880), Decimal(779299174080), Decimal("0.42284")),
    5: (Decimal(88366736882654697600), Decimal(7063261481362958280000), Decimal("0.38300411")),
    6: (Decimal(1080582891996915763252912268789760), Decimal(304351937651819767235871772704107520), Decimal("0.35910243")),
    7: (Decimal(134267527713152677457983522684235911142268448963584), Decimal(137077928808481869759712347602605837601908319657172480), Decimal("0.344418314818")),
}
out=[]
for n in range(2,8):
    dirname = "n2-strata-test" if n == 2 else f"n{n}-strata"
    d = json.loads((root / dirname / "aggregate.json").read_text())
    a1,a2,q=exact[n]
    A1=Decimal(str(d["A1_from_tau2"])); A2=Decimal(str(d["A2_sum"])); Q=Decimal(str(d["Q_e2_ratio"]))
    out.append({"n":n,"A1_ratio":str(A1/a1),"A2_ratio":str(A2/a2),"Q_error":str(Q-q),
                "min_A2_ess":d["min_A2_ess"],"max_A2_top_share":d["max_A2_top_share"]})
path=root/"small-grid-calibration.json";path.write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
