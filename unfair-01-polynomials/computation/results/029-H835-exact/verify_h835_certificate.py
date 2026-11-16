"""Independent top-offset 96-bit ball verifier; never imports the searcher."""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, struct, time
from verify_frontend import H, ALL, sha, verify
from reverse_ball_kernel import balls_from_boxes, top_coefficient

def req(ok, text):
    if not ok: raise ValueError(text)


def galois_check(path):
    # The local byte-identical integer module is independent of root/search code.
    from verify_modular_galois import verify as modular_verify
    raw = json.loads(Path(path).read_text()); rec = raw['certificates'][0]
    req(rec['coefficients_ascending'] == H and all(type(x) is int for x in rec['coefficients_ascending']),
        'fixed H8 galois input')
    req(H != list(reversed(H)), 'H8 nonreciprocal')
    out = modular_verify(rec)
    req(out['degree'] == 8 and out['status'] == 'full_symmetric_by_exact_modular_certificates', 'S8 modular status')
    return {'status': out['status'], 'galois_sha256': sha(path), 'nonreciprocal': True}


def raw_top(factors, indices, offset):
    lower, upper = [F(1)], [F(1)]
    for index in indices:
        factor = list(reversed(factors[index])); length = min(offset + 1, len(lower) + len(factor) - 1)
        lo, hi = [F(0)] * length, [F(0)] * length
        for i, (a, b) in enumerate(zip(lower, upper)):
            for j, (c, d) in enumerate(factor):
                if i + j >= length: break
                products = (a*c, a*d, b*c, b*d); lo[i+j] += min(products); hi[i+j] += max(products)
        lower, upper = lo, hi
    return lower[offset], upper[offset]


def validate_stream_metadata(meta, witness, canonical, balls_path, ledger_path):
    """Production validation for the compact u8 stream before any row is trusted."""
    count = len(canonical) // 4
    req(len(canonical) % 4 == 0, 'canonical u32 alignment')
    req(meta['format'] == 'H8-a3b5-top-offset-u8-v1' and meta['status'] == 'complete' and
        type(meta['bits']) is int and meta['bits'] == 128, 'meta format/status/bits')
    req(meta['balls_sha256'] == sha(balls_path) and meta['ledger_sha256'] == sha(ledger_path), 'input bindings')
    req(meta['canonical_count'] == count == 1742670 and meta['witness_bytes'] == len(witness) == 2 * count,
        'witness size')
    req(meta['witnesses_sha256'] == hashlib.sha256(witness).hexdigest(), 'witness hash')
    return count


def validate_row(word, side, offset, standard, factors, scale):
    """Production semantic check for one positional witness; standard returns None."""
    if word == standard:
        req((side, offset) == (255, 255), 'unique standard encoding')
        return None
    req(side in (0, 1) and 0 <= offset <= 64, 'u8 witness encoding')
    selected = word if side == 0 else word ^ ALL
    indices = [i for i in range(32) if selected >> i & 1]
    degree, (centre, radius) = top_coefficient(factors, indices, offset, scale)
    req(centre + radius < 0, 'strict negative top-offset ball')
    return indices, degree, centre, radius


def check(balls_path, ledger_path, certificate_path, witnesses_path, galois_path, record_path):
    balls = json.loads(Path(balls_path).read_text()); front, boxes, canonical = verify(balls, ledger_path)
    meta, witness = json.loads(Path(certificate_path).read_text()), Path(witnesses_path).read_bytes()
    validate_stream_metadata(meta, witness, canonical, balls_path, ledger_path)
    galois = galois_check(galois_path)
    scale, factors = balls_from_boxes(boxes, 96); standard = front['standard_mask']; negative = 0; minimum = None; samples=[]; seen=set(); start=time.monotonic()
    for pos, ((word,), pair) in enumerate(zip(struct.iter_unpack('<I', canonical), struct.iter_unpack('BB', witness)), 1):
        side, offset = pair
        validated = validate_row(word, side, offset, standard, factors, scale)
        if validated is None: continue
        indices, degree, centre, radius = validated
        margin = -(centre+radius); minimum = margin if minimum is None else min(minimum, margin); negative += 1
        if (offset,side) not in seen and len(samples)<12:
            seen.add((offset,side)); lo,hi=raw_top(boxes,indices,offset); req(hi<0,'raw Fraction sample'); samples.append({'position':pos-1,'offset':offset,'side':side,'actual_degree':degree,'raw_hi':str(hi)})
        if pos % 10000 == 0: record_path.with_suffix('.progress.json').write_text(json.dumps({'status':'verifying','checked':pos,'negative':negative})+'\n')
    req(negative == 1742669 and minimum is not None, 'negative/standard partition')
    result={'status':'all_H835_compact_witnesses_verified_by_independent_reversed_96bit_balls','negative_layouts':negative,'standard_layouts':1,'standard_mask':standard,'minimum_certified_negative_margin':str(F(minimum,scale)),'raw_fraction_sample_count':len(samples),'raw_fraction_samples':samples,'bits':96,'elapsed_seconds':time.monotonic()-start,'certificate_sha256':sha(certificate_path),'witnesses_sha256':sha(witnesses_path),'balls_sha256':sha(balls_path),'ledger_sha256':sha(ledger_path),'frontend_sha256':sha(Path(__file__).with_name('verify_frontend.py')),'reverse_ball_kernel_sha256':sha(Path(__file__).with_name('reverse_ball_kernel.py')),'helper_sha256':sha(Path(__file__).with_name('exact_complex.py')),'modular_galois_module_sha256':sha(Path(__file__).with_name('verify_modular_galois.py')),'verifier_sha256':sha(Path(__file__)),'galois':galois}
    record_path.write_text(json.dumps(result,indent=2)+'\n'); record_path.with_suffix('.progress.json').write_text(json.dumps({'status':'complete','checked':1742670,'negative':negative})+'\n'); print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for x in ('balls','ledger','certificate','witnesses','galois','record'): p.add_argument('--'+x,type=Path,required=True)
    a=p.parse_args(); check(a.balls,a.ledger,a.certificate,a.witnesses,a.galois,a.record)
