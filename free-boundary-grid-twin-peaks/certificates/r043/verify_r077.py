"""Recompute every n=7 sample multiplicity by a different exact algorithm."""
import concurrent.futures
import json
import sys
from collections import Counter
from pathlib import Path
from binary_cut_transfer import full_peaks, transfer

def count_chunk(path):
    data = json.loads(path.read_text())
    out = Counter()
    for word in data['samples']:
        x = transfer(7, word)
        assert x > 0
        peaks = full_peaks(7, word)
        assert len(peaks) == 2
        a, b = peaks
        assert abs(a//7-b//7)+abs(a%7-b%7)>7/3
        out[x] += 1
    return out

def main():
    source = Path(sys.argv[1])
    files = sorted((source/'raw').glob('*.json'))
    assert len(files) == 316
    histogram = Counter()
    with concurrent.futures.ProcessPoolExecutor(max_workers=12) as pool:
        for partial in pool.map(count_chunk, files, chunksize=4):
            histogram.update(partial)
    reference = json.loads((source/'summary.json').read_text())
    assert histogram == {int(k): v for k,v in reference['histogram_X'].items()}
    assert sum(histogram.values()) == 65536
    result = dict(n=7, chunks=len(files), samples=65536, exact_histogram_match=True,
                  recomputed_peaks_and_manhattan=True,
                  input_summary_sha256=__import__('hashlib').sha256((source/'summary.json').read_bytes()).hexdigest(),
                  histogram_X=dict(sorted(histogram.items())))
    Path(sys.argv[2]).write_text(json.dumps(result,indent=2)+'\n')
    print('65,536 words, aggregate histogram and peaks matched')

if __name__=='__main__':
    main()
