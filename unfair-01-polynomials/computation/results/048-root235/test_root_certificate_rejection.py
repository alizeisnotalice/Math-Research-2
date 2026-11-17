"""Regression: reject malformed coverage, disks, allocation, and negative index."""
import argparse,copy,json,subprocess,sys
from pathlib import Path

def main():
    p=argparse.ArgumentParser();p.add_argument('--certificate',type=Path,required=True);p.add_argument('--input',type=Path,required=True);p.add_argument('--diagnostics',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
    a.work.mkdir(parents=True,exist_ok=True);original=json.loads(a.certificate.read_text());script=Path(__file__).with_name('verify_root_disks_fixed.py');bad=a.work/'mutated-certificate.json';out=a.work/'unexpected-output.json'
    cases=['duplicate_task','zero_radius','wrong_root_center','bad_mask','constant_as_negative','optimized_python'];records=[]
    for label in cases:
        cert=copy.deepcopy(original)
        if label=='duplicate_task':cert['certificates'][0]['task']=19
        elif label=='zero_radius':cert['radius']['n']='0'
        elif label=='wrong_root_center':cert['certificates'][0]['centres_by_group'][0][0]['re']={'n':'100','d':'1'}
        elif label=='bad_mask':cert['certificates'][0]['mask_certificates'][0]['mask']=-1
        elif label=='constant_as_negative':cert['certificates'][0]['mask_certificates'][0]['negative_index']=0
        bad.write_text(json.dumps(cert));out.unlink(missing_ok=True)
        cmd=[sys.executable]+(['-O'] if label=='optimized_python' else [])+[str(script),'--certificate',str(bad),'--input',str(a.input),'--diagnostics',str(a.diagnostics),'--out',str(out)]
        r=subprocess.run(cmd,capture_output=True,text=True,timeout=180)
        if r.returncode==0 or out.exists():raise AssertionError('accepted invalid case '+label)
        records.append({'case':label,'returncode':r.returncode,'error_last_line':r.stderr.splitlines()[-1]})
    bad.unlink();(a.work/'rejection-tests.json').write_text(json.dumps({'status':'passed_all_six_rejections','cases':records},indent=2)+'\n');print('passed all six malformed/unsafe invocation cases')

if __name__=='__main__':main()
