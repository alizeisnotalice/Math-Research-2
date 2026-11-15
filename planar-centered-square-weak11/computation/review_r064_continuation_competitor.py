#!/usr/bin/env python3
"""Independent R64 accepted-chain/domain/endpoint audit and envelope checks."""
from fractions import Fraction as F
import hashlib,importlib.util,json,os
from pathlib import Path
import sys,time

HERE=Path(__file__).resolve().parent
ROOT=Path('/Volumes/data/CodexProjects/T2')
RUN=ROOT/'.local/compute/r064-continuation/run'
DRIVER=RUN.parent/'continue_certified_line.py'
PROOF=ROOT/'.local/workspaces/pcw-astra-covering/planar-centered-square-weak11/proofs/r064-overlapping-competitor-envelope.md'
MAIN_PROOF=PROOF.parent/'r064-source-state-mass-polytope-variation.md'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pin(p):return dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size)
def read(p):return json.loads(p.read_text())


def state_check(bb,mass,S,saved):
    eligible=[]
    for mask,b in sorted(bb.items()):
        m=sum(mass[i]for i in b['ids'])
        if m>b['d2']:eligible.append(dict(mask=mask,rank=mask.bit_count(),mass=str(m),span_squared=str(b['d2'])))
    assert eligible==saved['eligible_closures']
    general=max((F(e['mass'])for e in eligible if e['mask']!=S),default=F(0))
    overlap=max((F(e['mass'])for e in eligible if e['mask']!=S and e['mask']&S),default=F(0))
    assert str(general)==saved['eligible_competitor_max_mass_excluding_exact_S']
    assert str(overlap)==saved['competitor_overlap_max_mass']
    return overlap


def main():
    assert os.path.ismount('/Volumes/data')
    start=time.perf_counter()
    specmod=importlib.util.spec_from_file_location('independent_r64_geometry',HERE/'review_r064_mass_variation.py')
    audit=importlib.util.module_from_spec(specmod);specmod.loader.exec_module(audit)
    summary=read(RUN/'continuation-summary.json')
    assert summary['driver_sha256']==sha(DRIVER)=='5e2f9f7dd0a3c97402251c751eb13307af25fd9182a1bdbbda67236bc78c8ed6'
    assert summary['module_sha256']==sha(audit.PRODUCER)==audit.CODE
    assert summary['original_spec_sha256']==sha(audit.RUN/'input.json')==audit.INPUT
    initial=read(audit.RUN/'result.json');spec=read(audit.RUN/'input.json')
    assert summary['initial_result_sha256']==sha(audit.RUN/'result.json')
    assert summary['starting_objective']==initial['best_candidate']['objective']
    assert summary['original_seed_objective']==initial['seed_objective']
    points=[tuple(map(F,p))for p in spec['points']];bb=audit.boxes(points);S=spec['S_mask']
    mass=list(map(F,initial['best_candidate']['masses']));objective=summary['starting_objective']
    initial_overlap=state_check(bb,mass,S,summary['starting_state'])
    parent_hash=sha(audit.RUN/'result.json');records=[];artifacts=[pin(DRIVER),pin(RUN/'continuation-summary.json')]
    all_line_constraints=0;positive_saved_comparisons=0
    assert len(summary['steps'])==summary['budgets']['max_calls']==8
    assert summary['status']=='FINITE_CALL_BUDGET_EXHAUSTED'
    for i,step in enumerate(summary['steps'],1):
        ip,op=RUN/f'input-{i:02}.json',RUN/f'result-{i:02}.json'
        inp,result=read(ip),read(op)
        assert step['call']==i and step['accepted']and step['exit_status']==0
        assert step['input_sha256']==result['input_sha256']==sha(ip)
        assert step['parent_result_sha256']==inp['continuation_parent']['result_sha256']==parent_hash
        assert inp['continuation_parent']['original_input_sha256']==audit.INPUT
        assert step['output_sha256']==sha(op)
        assert result['source_sha256']==audit.CODE and result['status']=='CERTIFIED_INPUT_SPECIFIC_IMPROVEMENT'
        for k in ('points','S_indices','S_mask','mass_constraints','parameters','normalization','source_roles'):assert inp[k]==spec[k]
        assert list(map(F,inp['initial_masses']))==mass and result['seed_objective']==objective
        best=result['best_candidate'];newmass=list(map(F,best['masses']))
        check=audit.feasible(points,newmass,spec,bb)
        assert check['eligible_closures']==step['accepted_state']['exact_check']['strict_eligible_closures']
        assert check['nonterminal_margin']==step['accepted_state']['actual_nonterminal_margin']
        assert step['accepted_state']['connected_components']==[(1<<len(points))-1]
        assert step['accepted_state']['positive_source_count']==21
        overlap=state_check(bb,newmass,S,step['accepted_state'])
        assert best['objective']==step['accepted_objective']
        gain=F(best['objective']['lower'])-F(objective['upper'])
        assert gain>0 and str(gain)==step['strict_gain_lower']==best['certified_gain_lower']
        # Only accepted mass is independently checked here. Saved 8-proposal
        # scans are counted for scope, not geometrically recomputed.
        assert len(result['candidates'])==8
        assert best in result['candidates']
        for proposal in result['candidates']:
            proposal_mass=list(map(F,proposal['masses']))
            proposal_check=audit.feasible(points,proposal_mass,spec,bb)
            assert proposal['exact_check']['feasible'] and not proposal['exact_check']['errors']
            assert proposal_check['eligible_closures']==proposal['exact_check']['strict_eligible_closures']
            proposal_gain=F(proposal['objective']['lower'])-F(result['seed_objective']['upper'])
            assert str(proposal_gain)==proposal['certified_gain_lower']
            assert (proposal_gain>0)==proposal['certified_strict_improvement']
            all_line_constraints+=1;positive_saved_comparisons+=proposal_gain>0
        records.append(dict(call=i,parent_sha256=parent_hash,input_sha256=sha(ip),result_sha256=sha(op),
                            accepted_step=best['step'],constraint_check=check,eligible_catalog_exact_match=True,
                            overlap_competitor_mass=str(overlap),strict_gain_lower=str(gain),
                            gain_is_saved_interval_chain_arithmetic=True,geometry_recomputed_here=i==8))
        artifacts += [pin(ip),pin(op)]
        mass,objective,parent_hash=newmass,best['objective'],sha(op)
    assert mass==list(map(F,summary['final_masses']))and objective==summary['final_objective']
    final=audit.objective(points,mass,bb,S,72)
    for key,value in objective.items():
        if key!='scope':assert final[key]==value
    Lstar=F(records[-1]['overlap_competitor_mass'])
    assert Lstar==F(268181345,268435456)>F(3,4)
    final_gain_over_original_seed=F(final['lower'])-F(initial['seed_objective']['upper'])
    assert final_gain_over_original_seed>F(1,5)
    initial_positive=sum(bool(c['certified_strict_improvement'])for c in initial['candidates'])
    assert all_line_constraints==64 and initial_positive==8 and positive_saved_comparisons+initial_positive==42
    # Geometry-only envelope formula on a bounded exact rational grid.
    count=0
    for iu in range(16):
        for iv in range(16):
            for it in range(17):
                u,v,t=F(iu,16),F(iv,16),F(it,16)
                union=(u+t)*(v+t)+(1-u)*(1-v)-min(u+t,1-u)*min(v+t,1-v)
                assert 0<=union<=(1+t)**2
                # Both same-center rectangles fit the common square.
                assert max(u+t,1-u)<=1+t and max(v+t,1-v)<=1+t
                count+=1
    assert F(3)<F(49,16) # sqrt3<7/4, hence7/4+sqrt3<7/2.
    proof_hash=sha(PROOF)
    assert proof_hash=='3a9ae761af9c2f24fdfb49da3ecf291b4d2ba1d86cd98c651561580b6376ea28'
    assert sha(MAIN_PROOF)=='ef14cc7ee5b1384c241eb69bf78b6ae7639f1916afe58723e9187bef3a5cbdfe'
    out=dict(status='PASS_R64_CONTINUATION_ENDPOINT_AND_COMPETITOR_LEMMA',python=sys.version,executable=sys.executable,
             script_sha256=sha(Path(__file__)),elapsed_seconds=time.perf_counter()-start,source_pins=artifacts+[pin(PROOF),pin(MAIN_PROOF),pin(HERE/'review_r064_mass_variation.py')],
             accepted_chain=records,accepted_steps_checked=8,saved_intermediate_proposals=64,
             all_intermediate_proposal_geometry_recomputed=False,independent_geometry_cases_this_run=['continuation_endpoint'],
             continuation_line_mass_constraints_independently_checked=all_line_constraints,
             pilot_and_continuation_mass_feasible_points=all_line_constraints+len(initial['candidates']),
             pilot_and_continuation_positive_saved_interval_comparisons=positive_saved_comparisons+initial_positive,
             final_geometry=final,final_masses=list(map(str,mass)),final_L_star=str(Lstar),initial_best_L_star=str(initial_overlap),
             independently_certified_final_gain_over_original_seed=str(final_gain_over_original_seed),
             Lstar_is_conservative_eligible_closure_max_not_receiver_state_max=True,
             competitor_envelope_rational_cases=count,exact_corollary_comparison='3<49/16',
             stop_reason='Finite8-callbudget exhausted while final accepted gain positive; no local/global optimality or stationarity certificate.',
             scope='Independent all64 continuation proposal mass/domain constraints, all8 accepted-state catalog/hash/interval-chain audit, and endpoint2D area replay. No full geometry replay of all64proposal records; final accepted endpoint alone replayed here. No optimizer/MATLAB launch, C2 bound or3.5 conjecture proof.')
    target=HERE/'results/round064-continuation-competitor-independent.json';target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k not in ('source_pins','accepted_chain','final_masses')},indent=2))


if __name__=='__main__':main()
