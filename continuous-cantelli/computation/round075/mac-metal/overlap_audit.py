"""Recompute conservative four-device and native-command bounds for R075."""
from datetime import datetime, timedelta
import json
from pathlib import Path
import sys

BASE=Path(__file__).resolve().parent

def dt(t): return datetime.fromisoformat("2026-11-18T"+t+"Z")
def pairs(seq): return [x.isoformat() for x in seq]
def read(name): return json.loads((BASE/name).read_text())

def main():
    local_path=Path(sys.argv[1])
    metal=read('sync.json')
    off_min,off_max=-0.036053,0.219373  # MacBook clock minus local workstation.
    local_cpu=(dt('04:02:07.453'),dt('04:04:11.531'))
    local_metal=(dt('04:01:00.997'),dt('04:03:01.002'))
    remote_cpu=(dt('04:02:41.395347'),dt('04:04:41.004492'))
    safe_cpu=(remote_cpu[0]-timedelta(seconds=off_min),remote_cpu[1]-timedelta(seconds=off_max))
    remote_lifecycle=(dt('04:01:18'),dt('04:03:55'))
    start=max(local_cpu[0],safe_cpu[0],local_metal[0])
    finish=min(local_cpu[1],safe_cpu[1],local_metal[1])
    # The start marker precedes native initialization. Writing the result
    # follows the final sampled elapsed counter inside the result itself.
    ns=int((BASE/'sync.mtime.ns').read_text())
    mtime=datetime.fromtimestamp(ns/1e9,start.tzinfo)
    earliest=remote_lifecycle[0]-timedelta(seconds=off_max+0.001)
    latest=mtime-timedelta(seconds=metal['wallSeconds']+off_min-0.001)
    checkpoints=[]
    for ln in (BASE/'command-times.txt').read_text().splitlines():
        if ln.startswith('cycle='):
            d=dict(x.split('=',1) for x in ln.split())
            checkpoints.append((int(d['cycle']),float(d['wall']),float(d['gpu_seconds'])))
    lo=min(c for c in checkpoints if earliest+timedelta(seconds=c[1])>=start)
    hi=max(c for c in checkpoints if latest+timedelta(seconds=c[1])<=finish)
    mac_lower=hi[2]-lo[2]

    local=[]
    for ln in local_path.read_text().splitlines():
        if ln.startswith('utc=') or ln.startswith('end='):
            d=dict(x.split('=',1) for x in ln.split())
            t=datetime.fromisoformat(d.get('utc',d.get('end')))
            local.append((t,float(d['gpu_seconds'])))
    l0=max(t for t in local if t[0]<=start)
    l1=min(t for t in local if t[0]>=finish)
    local_idle=(l1[0]-l0[0]).total_seconds()-(l1[1]-l0[1])+0.002
    local_inner=[t for t in local if start<=t[0]<=finish]
    local_inner_lower=local_inner[-1][1]-local_inner[0][1]
    assert lo[0]==120 and hi[0]==144 and mac_lower>13.11 and local_idle<3.021
    summary={'rawMacCpuUTC':pairs(remote_cpu),'macCpuClockSafeUTC':pairs(safe_cpu),
      'macClockMinusLocalRangeSeconds':[off_min,off_max],
      'clockOffsetAssumedStable':True,
      'localCpuAll12Utc':pairs(local_cpu),'localMetalUtc':pairs(local_metal),
      'macMetalLifecycleUtc':pairs(remote_lifecycle),
      'fourWayTimeRangeUTC':pairs((start,finish)),
      'fourWaySeconds':(finish-start).total_seconds(),
      'macResultFileModificationTimeNs':ns,
      'macElapsedCounterZeroLocalBoundsUtc':pairs((earliest,latest)),
      'macContainedCycleCheckpoints':[lo,hi],
      'macContainedNativeCommandTimeSecondsLower':mac_lower,
      'localEnclosingCheckpoints':[[l0[0].isoformat(),l0[1]],[l1[0].isoformat(),l1[1]]],
      'localUnscheduledTimeSecondsUpper':local_idle,
      'localInnerCommandLower':local_inner_lower,
      'simultaneousGpuCommandSecondsViaEnclosingIdleLower':mac_lower-local_idle,
      'simultaneousGpuCommandSecondsViaBothContainedLower':mac_lower+local_inner_lower-(finish-start).total_seconds()}
    (BASE/'overlap.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(summary)

if __name__=='__main__': main()
