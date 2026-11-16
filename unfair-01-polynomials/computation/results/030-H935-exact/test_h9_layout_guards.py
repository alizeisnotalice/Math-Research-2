"""Run the actual independent layout checker against corrupted input files."""
import copy,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
def main():
    assert Path('/Volumes/data').is_mount()
    data=json.loads((HERE/'H9-real-necessary-layouts.json').read_text());results=[]
    for label in ['missing_layout','wrong_complement','wrong_root_selection']:
        bad=copy.deepcopy(data)
        if label=='missing_layout':bad['canonical_layouts'].pop()
        if label=='wrong_complement':bad['canonical_layouts'][5]['complement_mask_u37']^=1
        if label=='wrong_root_selection':bad['canonical_layouts'][5]['selected_group_indices'].append(36)
        with tempfile.TemporaryDirectory(prefix='h935-layout-',dir='/Volumes/data/CodexTemp/T2/unfair-01-polynomials') as tmp:
            t=Path(tmp);shutil.copyfile(HERE/'verify_h9_real_layout.py',t/'verify_h9_real_layout.py')
            (t/'H9-real-necessary-layouts.json').write_text(json.dumps(bad))
            run=subprocess.run([sys.executable,str(t/'verify_h9_real_layout.py')],capture_output=True,text=True)
            assert run.returncode!=0 and 'AssertionError' in run.stderr,(label,run.stdout,run.stderr)
            results.append({'test':label,'rejected_by_production_checker':True})
    (HERE/'H9-layout-guard-tests.json').write_text(json.dumps(results,indent=2)+'\n');print(results)
if __name__=='__main__':main()
