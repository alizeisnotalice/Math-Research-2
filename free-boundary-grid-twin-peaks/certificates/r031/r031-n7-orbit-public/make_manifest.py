#!/usr/bin/env python3
"""Write the public-package whitelist and SHA-256 manifest (excluding itself)."""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EXCLUDED = {"SHA256SUMS", "FILELIST.txt"}


def rel_files() -> list[Path]:
    return sorted((p for p in ROOT.rglob("*") if p.is_file() and p.relative_to(ROOT).as_posix() not in EXCLUDED),
                  key=lambda p: p.relative_to(ROOT).as_posix())


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main() -> None:
    files = rel_files()
    (ROOT / "FILELIST.txt").write_text("\n".join(p.relative_to(ROOT).as_posix() for p in files) + "\n")
    hashed = files + [ROOT / "FILELIST.txt"]
    (ROOT / "SHA256SUMS").write_text("".join(f"{digest(p)}  {p.relative_to(ROOT).as_posix()}\n" for p in hashed))
    print(f"wrote {len(files)} whitelisted files and {len(hashed)} SHA-256 entries")


if __name__ == "__main__":
    main()
