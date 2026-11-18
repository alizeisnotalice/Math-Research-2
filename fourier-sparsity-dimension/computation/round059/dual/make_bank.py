"""Deterministic affine partition bank for the R59 weighted covers."""
import argparse
import json
import struct
from pathlib import Path


class Generator:
    def __init__(self, initial):
        self.state=initial

    def get(self):
        self.state=(1664525*self.state+1013904223)&0xFFFFFFFF
        return self.state>>12


def bank():
    gen=Generator(0xA2805901)
    words=[]
    for n in (3,4):
        N=1<<n
        full=(1<<N)-1
        forms=[sum(((a&x).bit_count()%2)<<x for x in range(N)) for a in range(N)]
        for i in range(64):
            pending=[full]
            complete=[]
            while pending:
                a=pending.pop()
                depth=n-(a.bit_count().bit_length()-1)
                if a.bit_count()==1 or (depth and gen.get()%8<3):
                    complete.append(a)
                    continue
                choices=[row for row in forms[1:] if a & row and a & ~row]
                pick=choices[gen.get()%len(choices)]
                pending.extend([a&pick,a&~pick])
            complete.sort()
            assert sum(complete)==full
            words.extend([len(complete)]+complete+[0]*(16-len(complete)))
    return words


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    words=bank()
    args.output.write_bytes(struct.pack("<%dI"%len(words),*words))
