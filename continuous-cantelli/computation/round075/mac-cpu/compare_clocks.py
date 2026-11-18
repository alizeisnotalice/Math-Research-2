"""Bound relative post-run wall-clock offset with SSH round-trip receipts."""
import json
from pathlib import Path
import subprocess
import time


def main():
    checks=[]
    for _ in range(5):
        start=time.time_ns()
        p=subprocess.run(["ssh","-o","BatchMode=yes","-o","ConnectTimeout=5","macbook-math",
                          '~/.venvs/math/bin/python -c "import time; print(time.time_ns())"'],
                         text=True,capture_output=True,timeout=8)
        end=time.time_ns()
        if p.returncode==0:
            no=int(p.stdout.strip())
            checks.append(dict(local_before_ns=start,remote_ns=no,local_after_ns=end,
                               remote_minus_local_lower_ns=no-end,remote_minus_local_upper_ns=no-start))
    assert checks
    lower=max(x['remote_minus_local_lower_ns'] for x in checks)
    upper=min(x['remote_minus_local_upper_ns'] for x in checks)
    assert lower<=upper
    out=dict(post_run_remote_minus_local_ns=[lower,upper],checks=checks,
             model='offset bounds using each SSH request round-trip; measured after workload')
    Path("clock-comparison.json").write_text(json.dumps(out,indent=2)+"\n")
    print(lower,upper)


if __name__=='__main__':
    main()
