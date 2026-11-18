import json,sys,time,os
from pathlib import Path
from multiprocessing import get_context
from decimal import Decimal as D
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"negative"))
from certify_midpoint import certificate
def run(job):
    i,count,start,width,depth=job
    left=D('2.93595')+D(start+i)*D(width)
    right=left+D(width)
    t=time.time(); c=time.process_time()
    res=certificate(551,398,left,right,depth)
    return [i,os.getpid(),t,time.time(),time.process_time()-c,str(left),str(right),res['normalized_stripped_density']]
if __name__=='__main__':
    size,start=int(sys.argv[1]),int(sys.argv[2])
    with get_context('spawn').Pool(size) as pool:
        rows=pool.map(run,[(i,size,start,'0.000002',12) for i in range(size*2)],chunksize=1)
    Path(sys.argv[3]).write_text(json.dumps(rows,indent=2))
