#!/usr/bin/env python3
import json
from fractions import Fraction
from pathlib import Path

B=Path(__file__).parent
out={}
for n in (3,4):
    z=json.loads((B/f"n{n}.json").read_text())
    gain=sum((Fraction(int(v),n*n-int(k)-1) for k,v in z["rank_numerators"].items()),Fraction())
    out[str(n)]={"cover_records":z["cover_records"],"distinct_P":z["distinct_P"],
                 "pair_DPs":z["pair_DPs"],"cover_word_total":z["cover_word_total"],
                 "gain_num":str(gain.numerator),"gain_den":str(gain.denominator),
                 "average_num":str((gain/Fraction(int(z["cover_word_total"]))).numerator),
                 "average_den":str((gain/Fraction(int(z["cover_word_total"]))).denominator)}
assert out["3"]["gain_num"]=="51072" and out["3"]["gain_den"]=="1"
assert out["3"]["cover_word_total"]=="129024"
(B/"verified.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
