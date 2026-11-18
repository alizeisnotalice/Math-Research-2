"""Floating-point check of the piecewise symbolic primitive at the frozen five-body center."""
import json
from math import acos, atan2, cos, floor, hypot, pi, sin
from pathlib import Path

def evaluate(center):
    bodies = [(3,0,*center[:2]), (5,*center[2:5]), (7,*center[5:8]), (9,*center[8:])]
    end = 2*pi
    breaks = {0.,end}
    for n,angle,_,_ in bodies:
        for j in range(2*n):
            breaks.add((angle+(2*j+1)*pi/(2*n)) % end)
    def primitive(h,t):
        c,x,y = h
        return (c*c*t+2*c*(x*sin(t)-y*cos(t))+(x*x-y*y)*sin(2*t)/2-x*y*cos(2*t))/2
    total = 0.
    cells = 0
    bounds = sorted(breaks)
    for a,b in zip(bounds,bounds[1:]):
        mid = (a+b)/2
        rows = [(.5,0.,0.)]
        for n,angle,x,y in bodies:
            step = pi/n
            j = floor((mid-angle)/step+.5)
            sign = (-1)**j
            r = sign/(2*cos(step/2))
            rows.append(((1-sign)/2,x+r*cos(angle+j*step),y+r*sin(angle+j*step)))
        switches = {a,b}
        for i,(c,x,y) in enumerate(rows):
            for C,X,Y in rows[:i]:
                length = hypot(x-X,y-Y)
                if length and abs(C-c) <= length:
                    v = acos((C-c)/length)
                    for root in (atan2(y-Y,x-X)-v,atan2(y-Y,x-X)+v):
                        root %= end
                        if a < root < b:
                            switches.add(root)
        cut = sorted(switches)
        for l,r in zip(cut,cut[1:]):
            m = (l+r)/2
            h = max(rows,key=lambda p: p[0]+p[1]*cos(m)+p[2]*sin(m))
            total += primitive(h,r)-primitive(h,l)
            cells += 1
    return dict(cells=cells, area=total)

path = Path(__file__).with_name("five-body-contact-basin.json")
center = list(map(float,json.loads(path.read_bytes())["raw_center"]))
result = evaluate(center)
assert abs(result["area"]-.8372950062210252)<1e-10
print(json.dumps(result, indent=2))
