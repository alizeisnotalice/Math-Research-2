"""Integer checks for localization of mixed queries on matching selectors."""

import json
import random
import sys
from collections import Counter

SEED = 830021


def parity(value):
    return value.bit_count() % 2


def eliminate(value, pivots):
    mask = 0
    for bit in sorted(pivots, reverse=True):
        if value >> (bit - 1) & 1:
            row, combination = pivots[bit]
            value ^= row
            mask ^= combination
    return value, mask


def echelon(rows):
    pivots = {}
    for index, value in enumerate(rows):
        remainder, mask = eliminate(value, pivots)
        if remainder:
            pivots[remainder.bit_length()] = remainder, mask ^ (1 << index)
    return pivots


def independent_subset(rows):
    result = []
    for value in rows:
        if len(echelon(result + [value])) > len(result):
            result.append(value)
    return result


def rank(rows):
    return len(echelon(rows))


def affine_rank(points):
    return rank([2 * point + 1 for point in points])


def restrict(spectrum, rows, answer):
    pivots = echelon(rows)
    assert len(pivots) == len(rows)
    coefficients = Counter()
    for frequency, coefficient in spectrum.items():
        witness, mask = eliminate(frequency, pivots)
        coefficients[witness] += coefficient * (-1 if parity(mask & answer) else 1)
    return {frequency: coefficient for frequency, coefficient in coefficients.items()
            if coefficient}


def seed_walsh(dimension, function):
    size = 1 << dimension
    table = [function(point) for point in range(size)]
    return {frequency: value for frequency in range(size)
            if (value := sum(entry * (-1 if parity(point & frequency) else 1)
                             for point, entry in enumerate(table)))}


def bent(x, dimension):
    return -1 if sum(((x >> q) & 1)*((x >> (q+1)) & 1)
                     for q in range(0,dimension,2)) % 2 else 1


def edge_selector(n, edges, data_dim, seed):
    width = data_dim + 1
    spectrum = {}
    for i, (a,b) in enumerate(edges):
        offset = n + width*i
        for frequency in range(1 << n):
            address = (-1)**parity(frequency & a) + (-1)**parity(frequency & b)
            if not address:
                continue
            for k, coefficient in seed.items():
                spectrum[frequency | ((1 | (k << 1)) << offset)] = address*coefficient
    return spectrum


def block_support(word, n, width, blocks):
    word >>= n
    block_mask = (1 << width) - 1
    return tuple(i for i in range(blocks) if (word >> (i*width)) & block_mask)


def all_graph_values(mixed_rows, n, rhs):
    address_mask = (1 << n) - 1
    tuples = [(0,0)]
    for index, row in enumerate(mixed_rows):
        tuples += [(row ^ value, ((rhs >> index) & 1) ^ answer)
                   for value,answer in tuples]
    return {value & ~address_mask:(value, answer) for value,answer in tuples}


def two_endpoint_seed(seed, lower_rows, first_answer, second_answer):
    first = restrict(seed, lower_rows, first_answer)
    second = restrict(seed, lower_rows, second_answer)
    combined = {}
    for label in first.keys() | second.keys():
        a,b = first.get(label,0),second.get(label,0)
        if a+b:
            combined[2*label] = a+b
        if a-b:
            combined[2*label+1] = a-b
    return combined


def local_decomposition(n, edges, data_dim, seed, parent, mixed, answers, full, graph):
    width = data_dim + 1
    N = len(edges)
    address_mask = (1 << n) - 1
    local_rows = []
    local_answers = 0
    grouped_rows = [[] for _ in edges]
    for i in range(N):
        private_vectors = [private for private in graph if block_support(private,n,width,N) == (i,)]
        for private in independent_subset(private_vectors):
            row, answer = graph[private]
            assert not (row >> (n + i*width)) & 1
            if answer:
                local_answers |= 1 << len(local_rows)
            local_rows.append(row)
            grouped_rows[i].append((row,answer))
    intermediate = restrict(parent, local_rows, local_answers)
    selected = [[] for _ in edges]
    for frequency in intermediate:
        active = block_support(frequency, n, width, N)
        assert len(active) == 1
        selected[active[0]].append(frequency)

    r = affine_rank(seed)
    different_endpoints = 0
    for i,(a,b) in enumerate(edges):
        local = grouped_rows[i]
        local_seed_rows = [row >> (n + i*width+1) for row,_ in local]
        endpoint_answers = [
            sum((answer ^ parity((row & address_mask) & point)) << j
                for j,(row,answer) in enumerate(local))
            for point in (a,b)
        ]
        different_endpoints += endpoint_answers[0] != endpoint_answers[1]
        restricted_seed = two_endpoint_seed(seed, local_seed_rows, *endpoint_answers)
        assert len(selected[i]) == N*len(restricted_seed)
        projected_rank = rank([point >> n for point in selected[i]])
        assert 1 <= projected_rank <= min(r-len(local), affine_rank(restricted_seed))
        if not local:
            assert projected_rank == r

    local_pivots = echelon(local_rows)
    cross = independent_subset([eliminate(row,local_pivots)[0] for row in mixed])
    first_frequency = next(iter(intermediate))
    differences = [point ^ first_frequency for point in intermediate]
    omega = len(cross)+rank(differences)-rank(differences+cross)
    assert len(full) >= len(intermediate)-4*N*((1 << omega)-1)
    assert affine_rank(full) <= affine_rank(intermediate)-omega
    return {"local_queries": len(local_rows),
            "different_endpoint_answer_groups": different_endpoints,
            "effective_cross_queries": omega,
            "actual_cross_sparsity_loss": len(intermediate)-len(full)}


def check_query(n, edges, data_dim, seed, parent, a, mixed, answer):
    N = len(edges)
    width = data_dim + 1
    pure = [1 << q for q in range(a)]
    graph = all_graph_values(mixed,n,answer >> a)
    assert len(graph) == 1 << len(mixed)
    short = [x for x in graph if 1 <= len(block_support(x,n,width,N)) <= 2]
    marked = set(position for x in short for position in block_support(x,n,width,N))
    assert len(marked) <= 2*rank(short)

    belongs_to_cell = lambda point: point & ((1 << a)-1) == answer & ((1 << a)-1)
    occupancies = [int(belongs_to_cell(u))+int(belongs_to_cell(v)) for u,v in edges]
    f = sum(count==2 for count in occupancies)
    t = sum(count==1 for count in occupancies)
    M = 1 << (n-a)
    assert 2*f+t == M
    f_marked = sum(occupancies[i]==2 for i in marked)
    t_marked = sum(occupancies[i]==1 for i in marked)
    full = restrict(parent, pure+mixed, answer)
    assert len(full) >= len(seed)*(M//2*(f-f_marked)+M*(t-t_marked))

    stats = {"cases":1, "sum_marked_blocks":len(marked)}
    if not a:
        stats.update(local_decomposition(n,edges,data_dim,seed,parent,mixed,answer,full,graph))
    return stats, full


def random_edges(n, rng):
    points = list(range(1 << n))
    rng.shuffle(points)
    return list(zip(points[::2],points[1::2]))


def random_mixed(n, data_dim, blocks, w, rng):
    width = data_dim+1
    data_rows = []
    mixed = []
    while len(data_rows) < w:
        kind = rng.randrange(3)
        if kind == 0:
            i = rng.randrange(blocks)
            value = (2*rng.randrange(1,1 << data_dim)) << (i*width)
        elif kind == 1:
            i,j = rng.sample(range(blocks),2)
            value = ((1+2*rng.randrange(1 << data_dim)) << (i*width)) ^ (
                     (1+2*rng.randrange(1 << data_dim)) << (j*width))
        else:
            value = rng.randrange(1 << (width*blocks))
            if parity(value & sum(1 << (i*width) for i in range(blocks))):
                value ^= 1
        if rank(data_rows+[value]) > len(data_rows):
            data_rows.append(value)
            mixed.append((value << n) | rng.randrange(1 << n))
    return mixed


def main():
    rng = random.Random(SEED)
    result = {"python":sys.version.split()[0], "rng_seed":SEED, "audits":[]}
    for n, data_dim in ((2,2),(3,2),(4,2),(3,3)):
        edges = random_edges(n,rng)
        seed = seed_walsh(data_dim,
                         (lambda x:bent(x,data_dim)) if data_dim==2
                         else (lambda x:-1 if x==0 else 1))
        parent = edge_selector(n,edges,data_dim,seed)
        counts = Counter()
        for a in range(3):
            for w in range(1,4):
                for _ in range(8):
                    mixed = random_mixed(n,data_dim,len(edges),w,rng)
                    for answer in range(1 << (a+w)):
                        stats,_ = check_query(n,edges,data_dim,seed,parent,a,mixed,answer)
                        counts.update(stats)
        result["audits"].append({"address_dim":n,"seed_data_dim":data_dim,**counts})

    n,data_dim = 4,10
    edges = random_edges(n,rng)
    seed = seed_walsh(data_dim,lambda x:bent(x,data_dim))
    parent = edge_selector(n,edges,data_dim,seed)
    mixed = [(1 << n) ^ (1 << (n+data_dim+1))]
    parent_rank_plus_one = affine_rank(parent)+1
    rows = []
    for answer in (0,1):
        _,child = check_query(n,edges,data_dim,seed,parent,0,mixed,answer)
        child_rank_plus_one = affine_rank(child)+1
        assert len(seed) > 64*affine_rank(seed)
        assert len(parent) * child_rank_plus_one**2 < len(child) * parent_rank_plus_one**2
        rows.append({"support":len(child),"rank_plus_one":child_rank_plus_one})
    result["large_cross_only"] = {"parent_support":len(parent),
                                  "parent_rank_plus_one":parent_rank_plus_one,
                                  "children":rows}
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
