"""Exercise durable manifest ingestion and completion gating in an isolated root."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time

from scheduler import atomic, key_for, read_json, save

def main():
    root=Path(tempfile.mkdtemp(prefix="scheduler-smoke.",dir=sys.argv[1]))
    source=Path(__file__).with_name("scheduler.py")
    (root/"inbox").mkdir()
    save(root/"config.json",{"cpu_budget":int(sys.argv[2]),"metal_source":"unused"})
    def task(label,result):
        return dict(kind="cpu",stage="isolated-smoke",version="v1",scope=label,
                    input_sha256="1"*64,cpu=1,result={"path":"certificate.json","type":"json","equals":{"complete":True}},
                    argv=[sys.executable,"-B","-c",
                          "import pathlib,sys; pathlib.Path(sys.argv[1]).write_text(sys.argv[2])",
                          "{run}/certificate.json",json.dumps(result)])
    first=task("positive",{"complete":True})
    second=task("blocked-until-positive",{"complete":True})
    second["after"]=[key_for(first)]
    failure=task("wrong-result",{"complete":False})
    never=task("failed-parent-dependency",{"complete":True})
    never["after"]=[key_for(failure)]
    save(root/"inbox"/"batch.json",{"tasks":[second,first,failure,never]})
    save(root/"inbox"/"same-tasks.json",{"tasks":[second,first]})
    cmd=[sys.executable,"-B",str(source),"serve",str(root)]
    a=subprocess.Popen(cmd)
    try:
        deadline=time.monotonic()+45
        while time.monotonic()<deadline:
            hb=root/"heartbeat.json"
            if hb.exists() and read_json(hb)["pid"]==a.pid:
                break
            time.sleep(.25)
        else:
            raise RuntimeError("missing initial heartbeat")
        b=subprocess.run(cmd,capture_output=True,text=True,timeout=6)
        assert b.returncode==73,b.stdout+b.stderr
        while time.monotonic()<deadline:
            ledger=read_json(root/"ledger.json")
            states={row["task"]["scope"]:row["state"] for row in ledger["tasks"].values()}
            if states=={"positive":"complete","blocked-until-positive":"complete",
                        "wrong-result":"failed","failed-parent-dependency":"queued"}:
                break
            time.sleep(.25)
        else:
            raise RuntimeError(str(states))
        atomic(root/"STOP",b"smoke completion\n")
        assert a.wait(timeout=10)==0
        hb=read_json(root/"heartbeat.json")
        assert hb["stopping"] and hb["counts"]["complete"]==2
        print(json.dumps({"root":str(root),"singleton_exit":b.returncode,"states":states,"stop":"exit 0"}))
    finally:
        if a.poll() is None:
            a.terminate()
            a.wait(timeout=10)

if __name__=="__main__":
    main()
