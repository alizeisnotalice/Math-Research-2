"""Check fixed-angle cut containment and all exterior tangents using rational intervals."""
from fractions import Fraction as F
from pathlib import Path
from sigma_parameter import SIGMA_OUT
import contextlib
import hashlib
import io
import json
import runpy
import sys

base = Path(__file__).resolve().parent
sys.path.insert(0, str(base))


def run_script(name, *arguments):
    original = sys.argv
    path = base / name
    try:
        sys.argv = [str(path), *arguments]
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            namespace = runpy.run_path(str(path))
        return namespace, json.loads(output.getvalue())
    finally:
        sys.argv = original


sprague, sprague_result = run_script("certify_sprague_fixed.py")
fans, fans_result = run_script("certify_fans_in_p.py")
coordinates, coordinates_result = run_script("certify_ah_coordinates.py")
curvature, curvature_result = run_script("certify_ah_curvature.py")
eh, eh_result = run_script("audit_eh_geometry.py")

I, add, sub, mul, cross, dot, norm2, perp = (
    sprague[key] for key in ("I", "add", "sub", "mul", "cross", "dot", "norm2", "perp")
)
n, m = sprague["n"], sprague["m"]
half = I(F(1, 2))


def point(value):
    return tuple(I(item.lo, item.hi) for item in value)


def box(*points):
    return tuple(I(min(p[k].lo for p in points), max(p[k].hi for p in points))
                 for k in range(2))


def signs(vec, expected, label):
    for value, sign in zip(vec, expected):
        if sign:
            assert (value * sign).lo > 0, (label, value.out(), sign)
    return [v.out() for v in vec]


def out_points(value):
    return [item.out() for item in value]


def upper_cross(center_l, center_n):
    """The lower root has negative y whenever both centers are below zero."""
    q = sub(center_n, center_l)
    dsq = norm2(q)
    return sub(mul(half, add(center_l, center_n)),
               mul((1 - dsq / 4).sqrt() / dsq.sqrt(), perp(q)))


def support_slacks(value, exclude=()):
    supports = [(f"n{k}", u) for k, u in enumerate(n)]
    supports += [(f"m{k}", m[k]) for k in (3, 5)]
    result = {}
    for label, u in supports:
        if label in exclude:
            continue
        slack = half - dot(u, value)
        assert slack.lo > 0, (label, slack.out())
        result[label] = str(slack.lo)
    return result


base_ah = coordinates["d"]
C, G, E, Nminus, Y, Z, T_ah = (
    point(base_ah[k]) for k in ("C", "G", "E", "Nm", "Y", "Z", "T")
)
Tp = point(curvature["Tp"])
B1, F1, D1, E3, C3, F3, B3, K, Q, P, X, H, Ip = (
    sprague[k] for k in
    ("B1", "F1", "D1", "E3", "C3", "F3", "B3", "K", "Qp", "P", "X", "H", "Ip")
)
U, T, J = (point(eh[k]) for k in ("U", "T", "J"))
l_box, n_box = box(C, G), box(Nminus, E)
q = sub(n_box, l_box)
assert q[0].hi < 0 and l_box[1].hi < 0 and n_box[1].hi < 0
whole_r_box = upper_cross(l_box, n_box)
assert whole_r_box[1].lo > 0
anchor_disk_slacks = {"G": 1 - norm2(sub(Z, G)), "E3": 1 - norm2(sub(Z, E))}
assert all(slack.lo > 0 for slack in anchor_disk_slacks.values())
ah_halfplanes = support_slacks(whole_r_box, ("n1", "n2"))

# Locate the trimmed seam strictly within the Round 28 AS arcs.
seam = {}
for label, a, b, center in (
    ("Q--Y", Q, Y, E3), ("Y--X", Y, X, E3),
    ("X--Tp", X, Tp, G), ("Tp--P", Tp, P, G),
):
    u, v = sub(a, center), sub(b, center)
    c, d = cross(u, v), dot(u, v)
    assert c.lo > 0 and d.lo > 0, (label, c.out(), d.out())
    seam[label] = [str(c.lo), str(d.lo)]

# The EH arc coordinate extrema are at endpoints, since both radial
# vectors for each arc lie in quadrant III and span a minor arc.
for label, a, b, center in (
    ("U--J", U, J, C3), ("J--T", J, T, B3),
):
    u, v = sub(a, center), sub(b, center)
    signs(u, (-1, -1), label)
    signs(v, (-1, -1), label)
    assert cross(u, v).lo > 0 and dot(u, v).lo > 0
eh_box = box(U, J, T)
eh_halfplanes = support_slacks(eh_box, ("n3", "m3"))
eh_disk_slacks = {"U_to_B3": 1 - norm2(sub(U, B3)), "T_to_C3": 1 - norm2(sub(T, C3))}
assert all(slack.lo > 0 for slack in eh_disk_slacks.values())

fan_boxes = sprague["fan_boxes"]
cs, es_b, es_c, as_e, as_g = (fan_boxes[key] for key in ("CS", "ES-B", "ES-C", "AS-E", "AS-G"))
left = box(es_b, es_c, eh_box)
top = box(as_e, as_g, whole_r_box)
remote_margins = {
    "CS_right_of_ES_EH": cs[0].lo - left[0].hi,
    "AS_AH_above_CS": top[1].lo - cs[1].hi,
    "AS_AH_above_ES_EH": top[1].lo - left[1].hi,
}
assert all(value > 0 for value in remote_margins.values())

point_tangent = {
    "Q_start": perp(sub(Q, B1)),  # replaced below by the actual tangent
    "Q_circle": perp(sub(Q, E3)),
    "Y_from_E3": perp(sub(Y, E3)),
    "Y_to_C": perp(sub(Y, C)),
    "T_AH": perp(sub(T_ah, C)),
    "AH_segment": sub(Z, T_ah),
    "Tp_circle": perp(sub(Tp, G)),
    "P_circle": perp(sub(P, G)),
    "P_to_F1": sub(F1, P),
    "F1_to_H": sub(H, F1),
    "H_circle": perp(sub(H, C3)),
    "U_circle": perp(sub(U, C3)),
    "U_to_T": sub(T, U),
    "T_EH_circle": perp(sub(T, B3)),
    "I_circle": perp(sub(Ip, B3)),
    "I_to_E3": sub(E3, Ip),
    "E3_to_D1": sub(D1, E3),
    "D1_to_G": sub(G, D1),
    "G_circle": perp(sub(G, F3)),
    "K_circle": perp(sub(K, F3)),
    "K_to_C3": sub(C3, K),
    "C3_to_B1": sub(B1, C3),
}
point_tangent["Q_start"] = sub(Q, B1)
tangent_quadrants = {
    "Q_start": (-1, 1), "Q_circle": (-1, 1), "Y_from_E3": (-1, 1),
    "Y_to_C": (-1, -1), "T_AH": (-1, -1), "AH_segment": (-1, -1),
    "Tp_circle": (-1, -1), "P_circle": (-1, -1), "P_to_F1": (-1, -1),
    "F1_to_H": (0, -1), "H_circle": (0, -1), "U_circle": (1, -1),
    "U_to_T": (1, -1), "T_EH_circle": (1, -1), "I_circle": (1, -1),
    "I_to_E3": (1, -1), "E3_to_D1": (1, -1), "D1_to_G": (1, 1),
    "G_circle": (1, 1), "K_circle": (1, 1), "K_to_C3": (1, 1),
    "C3_to_B1": (0, 1),
}
tangent_signs = {key: signs(point_tangent[key], value, key)
                 for key, value in tangent_quadrants.items()}

join_cross, join_dot = {}, {}
for label, a, b in (
    ("Y", "Y_from_E3", "Y_to_C"), ("F1", "P_to_F1", "F1_to_H"),
    ("U", "U_circle", "U_to_T"), ("T_EH", "U_to_T", "T_EH_circle"),
    ("E3", "I_to_E3", "E3_to_D1"), ("D1", "E3_to_D1", "D1_to_G"),
    ("G", "D1_to_G", "G_circle"), ("C3", "K_to_C3", "C3_to_B1"),
    ("B1", "C3_to_B1", "Q_start"),
):
    a, b = point_tangent[a], point_tangent[b]
    c, d = cross(a, b), dot(a, b)
    assert c.lo > 0 and d.lo > 0, (label, c.out(), d.out())
    join_cross[label], join_dot[label] = str(c.lo), str(d.lo)

tangent_dot = {}
for label, a, b in (
    ("Q", "Q_start", "Q_circle"), ("T_AH", "T_AH", "AH_segment"),
    ("P", "P_circle", "P_to_F1"), ("H", "F1_to_H", "H_circle"),
    ("I", "I_circle", "I_to_E3"), ("K", "K_circle", "K_to_C3"),
):
    d = dot(point_tangent[a], point_tangent[b])
    assert d.lo > 0, (label, d.out())
    tangent_dot[label] = str(d.lo)

arc_tangents = {}
for a, b in (
    ("Q_circle", "Y_from_E3"), ("Y_to_C", "T_AH"),
    ("Tp_circle", "P_circle"), ("H_circle", "U_circle"),
    ("T_EH_circle", "I_circle"), ("G_circle", "K_circle"),
):
    c, d = cross(point_tangent[a], point_tangent[b]), dot(point_tangent[a], point_tangent[b])
    assert c.lo > 0 and d.lo > 0, (a, b, c.out(), d.out())
    arc_tangents[f"{a}--{b}"] = [str(c.lo), str(d.lo)]

result = {
    "passed": True,
    "sigma_degrees_interval": SIGMA_OUT,
    "python": sys.version,
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "dependencies_replayed": {
        "sprague": sprague_result["script_sha256"],
        "fans_in_p": fans_result["script_sha256"],
        "ah_coordinates": coordinates_result["script_sha256"],
        "ah_curvature": curvature_result["script_sha256"],
        "eh_geometry": eh_result["script_sha256"],
    },
    "AH_R_box": out_points(whole_r_box),
    "AH_anchor_Z_disk_slacks": {key: str(value.lo) for key, value in anchor_disk_slacks.items()},
    "AH_other_halfplane_slacks": ah_halfplanes,
    "AH_seam_cross_dot_lower_bounds": seam,
    "EH_box": out_points(eh_box),
    "EH_lens_endpoint_slacks": {key: str(value.lo) for key, value in eh_disk_slacks.items()},
    "EH_other_halfplane_slacks": eh_halfplanes,
    "remote_margins": {key: str(value) for key, value in remote_margins.items()},
    "join_cross_lower_bounds": join_cross,
    "join_dot_lower_bounds": join_dot,
    "tangent_forward_dot_lower_bounds": tangent_dot,
    "circle_ccw_short_arc_cross_dot_lower_bounds": arc_tangents,
    "tangent_component_intervals": tangent_signs,
}
print(json.dumps(result, indent=2))
