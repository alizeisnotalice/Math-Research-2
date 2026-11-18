"""Compare conditioned prefix ranks with the independent full n5 rank table."""
import json
from pathlib import Path
from merge_n11 import HERE, digest, full_square5_ranks, prefix_rank_counts, read, root_tables


def main():
    plan_file=HERE/'plan-n5-gate.json'
    result_file=HERE/'results-n5-gate/summary.json'
    plan=read(plan_file)
    book=read(result_file)
    assert digest(plan_file)==book['plan_sha256']
    ranks=prefix_rank_counts(plan,{(5,5,int(x['seen'])):int(x['completion']) for x in book['records']},root_tables())
    old=full_square5_ranks()
    checked=0
    for (w,h,u,c),row in ranks.items():
        first=abs(u%w-c)+(h-1-u//w)+1
        for j in set(range(1,5))|{first}:
            assert old[w,h,u,c][j]==row[j]
            checked+=1
    report={'completion_contexts':book['queried'],'bins_checked':checked,
            'n5_exact_independent_rank_archive_match':True,'query_ledger_sha256':digest(result_file)}
    (HERE/'n5-gate.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report)

if __name__=='__main__':main()
