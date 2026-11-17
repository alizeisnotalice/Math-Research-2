#!/usr/bin/env python3
"""Independent exact audit and ANF extraction for the R33 8-variable MM witness."""
import hashlib,json,platform,sys,time
from pathlib import Path

ROOT=Path(__file__).resolve().parent/'mm_completion_result.json'
SOURCE=Path(__file__).resolve().parents[1]/'round032/exceptional_signed_mask_result.json'
N=256
def fwht(a):
    a=list(a); h=1
    while h<len(a):
        for j in range(0,len(a),2*h):
            for i in range(j,j+h):
                x,y=a[i],a[i+h]; a[i],a[i+h]=x+y,x-y
        h*=2
    return a
def dot(x,y): return (x&y).bit_count()&1
def project(x,coordinates): return sum(((x>>c)&1)<<i for i,c in enumerate(coordinates))
def anf_coefficients(truth):
    a=truth[:]
    for bit in range(8):
        for mask in range(N):
            if mask&(1<<bit): a[mask]^=a[mask^(1<<bit)]
    return [i for i,c in enumerate(a) if c]
def monomial(mask,names):
    if mask==0:return '1'
    return '*'.join(names[i] for i in range(len(names)) if mask>>i&1)
def main():
    started=time.monotonic()
    root=json.loads(ROOT.read_text()); sol=root['solution']
    src=json.loads(SOURCE.read_text()); S=src['points']
    fhex=int(sol['f_hex'],16); ghex=int(sol['f_plus_h_hex'],16)
    f=[(fhex>>x)&1 for x in range(N)]; g=[(ghex>>x)&1 for x in range(N)]
    diff=[x for x in range(N) if f[x]^g[x]]
    assert diff==S and len(diff)==32 and all(x<128 for x in diff)
    def h_formula(x):
        b=[(x>>i)&1 for i in range(8)]
        return (1^b[7])*((b[0]*b[1]*b[3])^(b[0]*b[4]*b[5])^(b[1]*b[2]*b[6]))
    assert all((f[x]^g[x])==h_formula(x) for x in range(N))
    Wf=fwht([1-2*z for z in f]); Wg=fwht([1-2*z for z in g])
    assert set(abs(z) for z in Wf)=={16} and set(abs(z) for z in Wg)=={16}
    dual_switch=[i for i in range(N) if Wf[i]!=Wg[i]]
    dual_sums={a^b for a in dual_switch for b in dual_switch}
    dual_basis=[97,98,20];K={0}
    for a in dual_basis:K|={x^a for x in list(K)}
    assert len(dual_switch)==32 and len(dual_sums)==182 and len(K)==8 and K&dual_sums=={0}
    sums={a^b for a in diff for b in diff}
    assert len(sums)==128 and sums==set(range(128))
    data=tuple(sol['data_coords']); selector=tuple(sol['selector_coords']); pi=sol['pi']; rho=sol['rho']
    assert sorted(pi)==list(range(16))
    mm=[]
    for x in range(N):
        u=project(x,data); y=project(x,selector); mm.append(dot(u,pi[y])^rho[y])
    assert mm==f
    y_names=['y0','y1','y2','y3']; pi_anfs=[]
    for j in range(4):
        coeff=[(pi[y]>>j)&1 for y in range(16)]
        for bit in range(4):
            for mask in range(16):
                if mask&(1<<bit):coeff[mask]^=coeff[mask^(1<<bit)]
        pi_anfs.append([mask for mask,c in enumerate(coeff) if c])
    rv=rho[:]
    for bit in range(4):
        for mask in range(16):
            if mask&(1<<bit):rv[mask]^=rv[mask^(1<<bit)]
    f_anf=anf_coefficients(f); g_anf=anf_coefficients(g)
    out={'schema':'a28-r033-independent-mm-witness-audit-v1','scope':'independent exact verification of one supplied R33 witness; not a classification',
         'input_result_path':str(ROOT),'input_result_sha256':hashlib.sha256(ROOT.read_bytes()).hexdigest(),'switch_input_path':str(SOURCE),'switch_input_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
         'f_hex':sol['f_hex'],'f_plus_h_hex':sol['f_plus_h_hex'],'walsh_abs_f':sorted(set(abs(z) for z in Wf)),'walsh_abs_f_plus_h':sorted(set(abs(z) for z in Wg)),
         'dual_switch_support':dual_switch,'dual_switch_sumset_size':len(dual_sums),'dual_projection_kernel_basis':dual_basis,'switch_weight':len(diff),'switch_support':diff,'switch_sumset_size':len(sums),'switch_sumset_is_hyperplane_x7_zero':sums==set(range(128)),
         'mm_data_coordinates':list(data),'mm_selector_coordinates':list(selector),'pi':pi,'rho':rho,
         'pi_anf_masks':pi_anfs,'pi_anf_monomials':[[monomial(m,y_names) for m in a] for a in pi_anfs],
         'rho_anf_masks':[m for m,c in enumerate(rv) if c],'rho_anf_monomials':[monomial(m,y_names) for m,c in enumerate(rv) if c],
         'f_anf_masks':f_anf,'f_anf_monomials':[monomial(m,[f'x{i}' for i in range(8)]) for m in f_anf],'f_plus_h_anf_masks':g_anf,'f_plus_h_anf_monomials':[monomial(m,[f'x{i}' for i in range(8)]) for m in g_anf],
         'elapsed_seconds':time.monotonic()-started,'python':sys.version,'platform':platform.platform(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('independent-mm-witness-audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['walsh_abs_f','walsh_abs_f_plus_h','switch_weight','switch_sumset_size','mm_data_coordinates','mm_selector_coordinates','pi_anf_monomials','rho_anf_monomials','elapsed_seconds']},indent=2))
if __name__=='__main__': main()
