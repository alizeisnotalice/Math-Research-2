#!/usr/bin/env python3
"""Read-only research inventory; hashes identify copies, never proof validity."""
import argparse
import csv
import hashlib
import json
import os
from collections import Counter, defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('/Volumes/data/CodexProjects/T'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    if not source.is_dir() or output == source or source in output.parents:
        parser.error('Source must exist; output must be outside the legacy tree.')
    output.mkdir(parents=True, exist_ok=True)
    roots = sorted(p for p in source.iterdir() if p.is_dir() and not p.is_symlink()
                   and (p.name.startswith(('Math-Research', 'cube-', 'centered-cube', 'planar-centered'))
                        or p.name == '八维中心方体论文_20261027'))
    skipped = {'.git', '__pycache__', '.cache', 'node_modules', '.venv', 'venv', 'build'}
    text_ext = {'.md', '.tex', '.bib', '.py', '.m', '.wl', '.wls', '.json', '.csv'}
    digests = defaultdict(list)
    rows, errors = [], []
    with (output / 'files.csv').open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(['path_relative_to_T', 'bytes', 'sha256_text_up_to_2MiB'])
        for root in roots:
            counts, total, text_count = Counter(), 0, 0
            for directory, dirs, files in os.walk(root, followlinks=False):
                dirs[:] = sorted(d for d in dirs if d not in skipped and not (Path(directory)/d).is_symlink())
                for name in sorted(files):
                    p = Path(directory) / name
                    if p.is_symlink():
                        continue
                    try:
                        size = p.stat().st_size
                        total += size
                        counts[p.suffix.lower() or '(none)'] += 1
                        digest = ''
                        if p.suffix.lower() in text_ext and size <= 2 * 1024**2:
                            digest = hashlib.sha256(p.read_bytes()).hexdigest()
                            digests[digest].append(str(p.relative_to(source)))
                            text_count += 1
                        writer.writerow([p.relative_to(source), size, digest])
                    except OSError as exc:
                        errors.append({'path': str(p.relative_to(source)), 'error': type(exc).__name__})
            rows.append({'directory': root.name, 'files': sum(counts.values()), 'bytes': total,
                         'hashed_text_files': text_count, 'extensions': dict(sorted(counts.items()))})
    groups = [v for _, v in sorted(digests.items()) if len(v) > 1]
    summary = {'roots': len(rows), 'files': sum(r['files'] for r in rows),
               'hashed_text_files': sum(r['hashed_text_files'] for r in rows),
               'unique_text_hashes': len(digests), 'duplicate_text_groups': len(groups),
               'duplicate_text_copies': sum(len(g)-1 for g in groups),
               'excluded_directory_names': sorted(skipped), 'errors': errors,
               'note': 'Indexing and byte equality only; not semantic reading or mathematical verification.'}
    for name, data in [('summary.json', summary), ('roots.json', rows), ('duplicates.json', groups)]:
        (output/name).write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == '__main__':
    main()
