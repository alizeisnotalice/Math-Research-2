"""Check transferred range coverage and numerical metadata for the Metal run."""

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path


def digest(path):
    with Path(path).open("rb") as payload:
        return sha256_stream(payload)


def sha256_stream(payload):
    hasher=sha256()
    for block in iter(lambda:payload.read(1024*1024),b""):
        hasher.update(block)
    return hasher.hexdigest()


def main():
    assert __debug__
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs",type=Path)
    parser.add_argument("final",type=Path)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    packet=json.loads(args.inputs.read_bytes())
    final=json.loads(args.final.read_bytes())
    assert final["status"]=="complete_numerical_coverage"
    manifest=Path(final["handoff"])
    assert digest(manifest)==final["sha256"]
    doc=json.loads(manifest.read_bytes())
    assert doc["source_inputs_sha256"]==digest(args.inputs)
    by_task={task["id"]:task["masks"] for task in packet["tasks"]}
    current={task:0 for task in by_task}
    count=Counter()
    small=None
    for part in doc["layers"]:
        receipt=Path(part["receipt_path"])
        assert digest(receipt)==part["receipt_sha256"]
        evidence=json.loads(receipt.read_bytes())
        assert evidence["status"]=="complete" and evidence["coverage_credit"] is True
        task=evidence["task"]
        lo,hi=evidence["requested_range"]
        assert lo==current[task] and lo<hi<=by_task[task]
        if evidence["kind"]=="recovered-degree54-stream-receipt":
            assert len(evidence["remote_final_sha256"])==64
            count["recovered_segments"]+=1
        else:
            assert evidence["kind"]=="single_process_metal_stream_receipt"
            assert evidence["verified_prefix_end"]==hi
        stream=Path(part["path"])
        assert digest(stream)==part["sha256"]==evidence["ledger_sha256"]
        batches=0
        with stream.open("rb") as payload:
            for line in payload:
                event=json.loads(line)
                start,end=event["range"]
                assert event["task"]==task and event["coverage_credit"] is True
                assert start==current[task] and end==min(hi,start+131072)
                assert event["finite_count"]==event["count"]==end-start
                candidates=event["candidates_cpu_reconstructed"]
                assert len(candidates)==4 and len({x["mask"] for x in candidates})==4
                for x in candidates:
                    assert start<=x["mask"]<end
                    assert isfinite(x["cpu_direct_score"]) and isfinite(x["gpu_float_score"])
                    assert isfinite(x["target_product_max_abs_error"])
                    count["numeric_unknown"]+=int(x["near_threshold_numeric_unknown"])
                    if small is None or x["cpu_direct_score"]<small["cpu_direct_score"]:
                        small={"task":task,"mask":x["mask"],"cpu_direct_score":x["cpu_direct_score"],
                               "product_residual":x["target_product_max_abs_error"]}
                count["retained_samples"]+=len(candidates)
                count["batches"]+=1
                count["mask_count"]+=end-start
                batches+=1
                current[task]=end
        assert current[task]==hi and batches==evidence["verified_batches"]
    assert current==by_task and list(doc["task_next"].values())==[by_task[int(task)] for task in doc["task_next"]]
    assert count["mask_count"]==sum(by_task.values())
    result={"status":"numeric_coverage_verified","inputs_sha256":digest(args.inputs),"handoff_sha256":digest(manifest),
            "task_count":len(by_task),"segments":len(doc["layers"]),"counts":dict(count),"minimum_cpu_score":small,
            "scope":"Range and saved metadata audit only; scores were not recomputed at all masks."}
    with args.out.open("x") as output:
        json.dump(result,output,indent=2)
        output.write("\n")
    print(json.dumps(result))


if __name__=="__main__":
    main()
