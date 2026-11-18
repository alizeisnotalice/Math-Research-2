"""Integer spectrum and rank checks for arbitrary signed selector masks."""

import json
import random
import sys
from collections import Counter
from fractions import Fraction

SEED = 840021


def walsh(values):
    coefficients = list(values)
    width = 1
    while width < len(coefficients):
        for start in range(0, len(coefficients), 2*width):
            for i in range(start, start+width):
                j = i+width
                left,right = coefficients[i],coefficients[j]
                coefficients[i],coefficients[j] = left+right,left-right
        width *= 2
    return {point: value for point,value in enumerate(coefficients) if value}


def independent(vectors):
    chosen = []
    pivots = {}
    for vector in vectors:
        remainder = vector
        while remainder:
            bit = remainder.bit_length()
            if bit not in pivots:
                pivots[bit] = remainder
                chosen.append(vector)
                break
            remainder ^= pivots[bit]
    return chosen


def rank(vectors):
    return len(independent(vectors))


def span(basis):
    points = [0]
    for vector in basis:
        points += [point ^ vector for point in points]
    return points


def affine_rank(frequencies):
    return rank([2*frequency+1 for frequency in frequencies])


def phase_restriction(cell, phase):
    anchor = min(cell)
    names = span(independent([p^anchor for p in cell]))
    assert {p^anchor for p in names} == set(cell)
    return [phase[p^anchor] for p in names]


def mask_spectrum(cell, phase):
    members = set(cell)
    return walsh([sign if position in members else 0
                  for position,sign in enumerate(phase)])


def private_data(address_dim, cells, phase, children):
    child_dim = len(children[0]).bit_length()-1
    address_differences = []
    width = child_dim+1
    parent = {}
    mask_width = []
    kid_width = []
    kid_rank = []
    for i,(cell,child) in enumerate(zip(cells,children)):
        address = mask_spectrum(cell,phase)
        data = walsh(child)
        mask_width.append(len(address))
        kid_width.append(len(data))
        kid_rank.append(affine_rank(data))
        anchor = next(iter(address))
        address_differences += [p^anchor for p in address]
        for a,ac in address.items():
            for beta,bc in data.items():
                guard_label = (2*beta+1) << (address_dim+i*width)
                parent[a | guard_label] = ac*bc
    return parent,mask_width,kid_width,kid_rank,rank(address_differences)


def private_truth(address_dim, cells, phase, children):
    child_dim = len(children[0]).bit_length()-1
    private_width = child_dim+1
    address_mask = (1 << address_dim)-1
    child_mask = (1 << child_dim)-1
    labels = {p:i for i,cell in enumerate(cells) for p in cell}
    for physical in range(1 << (address_dim+private_width*len(cells))):
        u = physical & address_mask
        i = labels[u]
        block = physical >> (address_dim+i*private_width)
        guard = (-1)**(block&1)
        yield phase[u]*guard*children[i][(block >> 1)&child_mask]


def check_private(address_dim, cells, phase, children, check_truth):
    spectrum,widths,sizes,ranks,t = private_data(address_dim,cells,phase,children)
    parent_rank = affine_rank(spectrum)
    assert len(spectrum) == sum(k*s for k,s in zip(widths,sizes))
    assert parent_rank == sum(ranks)+t
    if check_truth:
        kid_dim = len(children[0]).bit_length()-1
        common_factor = 2**((len(cells)-1)*(kid_dim+1)+1)
        actual = walsh(private_truth(address_dim,cells,phase,children))
        assert actual == {frequency:common_factor*value for frequency,value in spectrum.items()}
    gamma = sum(Fraction(1,k) for k in widths)
    assert gamma <= 1
    extra = t+1-len(cells)
    if gamma == 1:
        assert extra <= 0
    norm = sum(Fraction((r+1)**2,k*s) for k,s,r in zip(widths,sizes,ranks))
    assert (parent_rank+1-max(extra,0))**2 <= len(spectrum)*norm
    return parent_rank, len(spectrum), t, gamma


def weighted_affine_bound(cells, phase, children):
    total = Fraction(0)
    for cell,child in zip(cells,children):
        intrinsic = walsh(phase_restriction(cell,phase))
        data = walsh(child)
        r_phase,r_data = affine_rank(intrinsic),affine_rank(data)
        product_cost = Fraction((r_phase+r_data)**2,len(intrinsic)*len(data))
        assert product_cost <= max(Fraction((r_phase+1)**2,len(intrinsic)),
                                   Fraction((r_data+1)**2,len(data)))
        total += Fraction(len(cell),len(phase))*product_cost
    return total


def random_partition(address_count, cells, rng):
    order = list(range(address_count))
    rng.shuffle(order)
    strips = [[p] for p in order[:cells]]
    for point in order[cells:]:
        strips[rng.randrange(cells)].append(point)
    return strips


def affine_refining_pairs(cells, address_count):
    labels = {p:i for i,cell in enumerate(cells) for p in cell}
    pieces = []
    parents = []
    for a in range(0,address_count,2):
        if labels[a] == labels[a+1]:
            pieces.append([a,a+1])
            parents.append(labels[a])
        else:
            pieces.extend([[a],[a+1]])
            parents.extend([labels[a],labels[a+1]])
    return pieces,parents


def check_identification(cells, phase, children):
    address_count = len(phase)
    address_dim = address_count.bit_length()-1
    pieces,parents = affine_refining_pairs(cells,address_count)
    expanded_children = [children[i] for i in parents]
    expanded_rank,expanded_width,_,_ = check_private(
        address_dim,pieces,phase,expanded_children,False
    )
    summed_coefficients = Counter()
    multiplicity = Counter()
    for j,cell in enumerate(pieces):
        for alpha,ac in mask_spectrum(cell,phase).items():
            for beta,bc in walsh(expanded_children[j]).items():
                frequency = beta*address_count+alpha
                summed_coefficients[frequency] += ac*bc
                multiplicity[frequency] += 1
    final = {point:value for point,value in summed_coefficients.items() if value}
    labels = {p:i for i,cell in enumerate(cells) for p in cell}
    final_truth = [phase[p]*children[labels[p]][data]
                   for data in range(len(children[0])) for p in range(address_count)]
    assert final == walsh(final_truth)
    assert affine_rank(final) <= expanded_rank
    repetitions = sum(n-1 for n in multiplicity.values())
    cancellations = sum(point not in final for point in multiplicity)
    assert expanded_width-len(final) == repetitions+cancellations
    penalty = Fraction(expanded_width,len(final))
    bound = weighted_affine_bound(pieces,phase,expanded_children)
    assert Fraction((affine_rank(final)+1)**2,len(final)) <= penalty*bound
    return repetitions,cancellations


def check_phase_translation(address_dim, child, rng):
    address_count = 1 << address_dim
    signs = [rng.choice((-1, 1)) for _ in range(address_count)]
    translations = [rng.randrange(len(child)) for _ in range(address_count)]
    actual = walsh([signs[u]*child[data ^ translations[u]]
                    for data in range(len(child)) for u in range(address_count)])
    predicted = {}
    differences = []
    kid_spectrum = walsh(child)
    noncharacter_count = 0
    for beta,value in kid_spectrum.items():
        mask = walsh([sign * (-1)**((beta & shift).bit_count() % 2)
                      for sign,shift in zip(signs,translations)])
        noncharacter_count += len(mask)>1
        anchor = next(iter(mask))
        differences.extend(point^anchor for point in mask)
        for alpha,coefficient in mask.items():
            predicted[beta*address_count+alpha] = value*coefficient
    assert actual == predicted
    assert affine_rank(actual) == affine_rank(kid_spectrum)+rank(differences)
    if len(kid_spectrum)==len(child) and noncharacter_count:
        assert 2*noncharacter_count >= len(child)
        assert 2*len(actual) >= 5*len(child)


def main():
    rng = random.Random(SEED)
    catalogue = [
        [1,1,1,1], [1,-1,1,-1], [1,1,1,-1], [-1,1,1,1],
    ]
    totals = Counter()
    for index in range(192):
        cells = random_partition(8, 2+index%2, rng)
        phase = [rng.choice((-1,1)) for _ in range(8)]
        children = [rng.choice(catalogue) for _ in cells]
        r,s,_,gamma = check_private(3,cells,phase,children,True)
        repetitions,cancellations = check_identification(cells,phase,children)
        totals["arbitrary_signed_masks"] += 1
        totals["uncertainty_equality"] += gamma == 1
        totals["identified_repetitions"] += repetitions
        totals["cancelled_fibres"] += cancellations
        totals["tested_actual_spectral_points"] += s

        affine_cells = [list(range(4)),list(range(4,8))]
        affine_children = children[:2]
        r,s,_,_ = check_private(3,affine_cells,phase,affine_children,True)
        assert Fraction((r+1)**2,s) <= weighted_affine_bound(
            affine_cells,phase,affine_children
        )
        totals["affine_signed_masks"] += 1
    for index in range(128):
        check_phase_translation(3+index%2,rng.choice(catalogue),rng)
        totals["arbitrary_phase_translations"] += 1

    line_a,plane,line_b = [0,4],[1,3,5,7],[2,6]
    phase = [(-1)**((value >> 2)&1) if value&1 else 1
             for value in range(8)]
    parts = [mask_spectrum(cell,phase) for cell in (line_a,plane,line_b)]
    assert [len(part) for part in parts] == [4,2,4]
    assert not (parts[0].keys() & parts[1].keys())
    assert len(mask_spectrum(line_a+plane,phase)) == 6
    totals["nonaffine_zero_penalty_template"] = 1
    print(json.dumps({"python":sys.version.split()[0],"seed":SEED,**totals},indent=2))


if __name__ == "__main__":
    main()
