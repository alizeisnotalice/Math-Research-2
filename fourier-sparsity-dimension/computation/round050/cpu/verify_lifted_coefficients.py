"""Independently verify the three explicit spectra and an affine change of variables."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from affine_incidence_audit import ADDR, AMBIENT, Partition, fwht, parity, rank


def pullback(points, basis):
    coordinates = [0]
    for b in basis:
        coordinates.extend(a ^ b for a in coordinates[:])
    assert len(set(coordinates)) == len(coordinates)
    return coordinates[points]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_cases", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    data = json.loads(args.source_cases.read_bytes())
    dim = AMBIENT.bit_length() - 1
    bases = [1 << i for i in range(dim)]
    for i, j in ((0, 3), (3, 2), (1, 4), (5, 0), (6, 2), (2, 4)):
        bases[i] ^= bases[j]
    assert rank(bases) == dim
    affine_offset = 53
    physical_map = [pullback(x, bases) ^ affine_offset for x in range(AMBIENT)]
    inverse_map = {image: x for x, image in enumerate(physical_map)}
    result = {}
    for name, case in data["examples"].items():
        geo = Partition(case["blocks"])
        truth_mask = int(case["truth_mask"], 16)
        ordinary = fwht(1 - 2 * ((truth_mask >> x) & 1) for x in range(AMBIENT))
        lifted_coefficients = [0] * AMBIENT
        term_supports = []
        for i, shape in enumerate(case["leaf_phases_q_r_linear_constant"]):
            _, _, _, sign = geo.make_leaf(i, shape)
            addresses, address_basis = geo.coordinates[i]
            small = fwht(sign)
            terms = []
            for freq in range(AMBIENT):
                local = geo.project(freq, address_basis, ADDR.bit_length() - 1)
                value = (1 - 2 * parity(freq & addresses[0])) * small[local]
                lifted_coefficients[freq] += value
                terms.append(value)
            term_supports.append(terms)
        assert ordinary == lifted_coefficients

        rotated_truth = fwht(1 - 2 * ((truth_mask >> z) & 1) for z in physical_map)
        expected_rotated = [0] * AMBIENT
        for old_freq, val in enumerate(ordinary):
            to = sum(parity(old_freq & base) << i for i, base in enumerate(bases))
            expected_rotated[to] = (1 - 2 * parity(old_freq & affine_offset)) * val
        assert rotated_truth == expected_rotated
        common = geo.address_basis + [1 << i for i in range(3, dim)]
        rotated_common = [inverse_map[v] ^ inverse_map[0] for v in common]
        assert rank(rotated_common) == len(common)
        counts = Counter(
            sum(parity(v & f) << i for i, v in enumerate(rotated_common))
            for f, value in enumerate(rotated_truth) if value
        )
        assert dict(counts) == {int(beta): value["band_support"] for beta, value in case["band_details"].items()}
        verified_private = 0
        for freq, value in enumerate(ordinary):
            beta = str(geo.projection[freq])
            if beta not in case["band_details"]:
                continue
            participants = case["band_details"][beta]["participants_local_preimage_counts"]
            if len(participants) == 1:
                i, _ = participants[0]
                assert value == term_supports[i][freq]
                assert all(t[freq] == 0 for j, t in enumerate(term_supports) if i != j)
                verified_private += bool(value)
        result[name] = {"all_lifted_coefficients_matched": AMBIENT, "private_lifted_coefficients": verified_private,
                        "rotated_L_basis": rotated_common, "beta_sparsity_after_change": dict(counts)}

    out = {
        "input_sha256": hashlib.sha256(args.source_cases.read_bytes()).hexdigest(),
        "basis_columns_of_affine_change": bases, "affine_offset": affine_offset, "examples": result,
    }
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
