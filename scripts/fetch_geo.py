#!/usr/bin/env python3
"""Download GEO series supplementary files and series-matrix metadata.

    python scripts/fetch_geo.py GSE172258 GSE81700 ...      # explicit list
    python scripts/fetch_geo.py --all                         # every dataset in config/datasets.yaml
    python scripts/fetch_geo.py --all --list-only             # only write listings (fast)

For each series it writes, under data/raw/<GSE>/:
    listing.tsv        every supplementary/matrix file GEO offers, with size
    manifest.tsv       what was downloaded, with SHA-256 and UTC date
    skipped.tsv        files above --max-gb (re-run with a larger cap if needed)
and unpacks *_RAW.tar archives into RAW/.

Standard library only; respects HTTPS_PROXY. Re-running skips files already complete.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import sys
import tarfile
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series"
ROOT = Path(__file__).resolve().parents[1]
HREF = re.compile(r'<a href="([^"?/][^"]*)">[^<]*</a>\s*(\S+\s+\S+)?\s*(\S+)?')
UNITS = {"K": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}


def stem(gse: str) -> str:
    digits = gse[3:]
    return f"GSE{digits[:-3]}nnn" if len(digits) > 3 else "GSEnnn"


def get(url: str, retries: int = 3) -> bytes:
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except Exception:
            if attempt == retries:
                raise
            time.sleep(3 * attempt)
    raise RuntimeError("unreachable")


def parse_size(text: str | None) -> float:
    if not text:
        return float("nan")
    m = re.fullmatch(r"([\d.]+)([KMGT]?)", text.strip())
    if not m:
        return float("nan")
    return float(m.group(1)) * UNITS.get(m.group(2), 1)


def listing(gse: str) -> list[tuple[str, str, float]]:
    """(subdir, filename, approx_bytes) for suppl/ and matrix/."""
    rows = []
    for sub in ("suppl", "matrix"):
        url = f"{GEO}/{stem(gse)}/{gse}/{sub}/"
        try:
            html = get(url).decode("utf-8", "replace")
        except Exception as exc:
            print(f"  [warn] {url}: {exc}")
            continue
        for m in HREF.finditer(html):
            name = m.group(1)
            if name.startswith(("http", "#")) or name.endswith("/"):
                continue
            rows.append((sub, name, parse_size(m.group(3))))
    return rows


def download(url: str, dest: Path) -> None:
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, tmp.open("wb") as fh:
                shutil.copyfileobj(r, fh, length=1 << 20)
            tmp.replace(dest)
            return
        except Exception as exc:
            print(f"  [retry {attempt}] {dest.name}: {exc}")
            time.sleep(5 * attempt)
    raise RuntimeError(f"failed to download {url}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(gse: str, out_root: Path, max_bytes: float, list_only: bool) -> None:
    print(f"== {gse}")
    out = out_root / gse
    out.mkdir(parents=True, exist_ok=True)
    files = listing(gse)
    with (out / "listing.tsv").open("w") as fh:
        fh.write("subdir\tfile\tapprox_bytes\n")
        fh.writelines(f"{s}\t{n}\t{b:.0f}\n" for s, n, b in files)
    print(f"  {len(files)} file(s) listed")
    if list_only:
        return

    manifest = out / "manifest.tsv"
    if not manifest.exists():
        manifest.write_text("file\turl\tbytes\tsha256\tdownloaded_utc\n")
    skipped = []
    for sub, name, size in files:
        dest = out / name
        url = f"{GEO}/{stem(gse)}/{gse}/{sub}/{name}"
        if size == size and size > max_bytes:  # size == size: not NaN
            skipped.append((name, size))
            print(f"  [skip] {name} (~{size / 1e9:.1f} GB > cap)")
            continue
        if dest.exists() and dest.stat().st_size > 0:
            print(f"  [have] {name}")
            continue
        print(f"  [get ] {name} (~{size / 1e6:.0f} MB)")
        download(url, dest)
        with manifest.open("a") as fh:
            stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            fh.write(f"{name}\t{url}\t{dest.stat().st_size}\t{sha256(dest)}\t{stamp}\n")
        if name.endswith("_RAW.tar"):
            raw = out / "RAW"
            raw.mkdir(exist_ok=True)
            with tarfile.open(dest) as tar:
                tar.extractall(raw, filter="data")
            print(f"  [untar] → {raw}")
    if skipped:
        with (out / "skipped.tsv").open("w") as fh:
            fh.write("file\tapprox_bytes\n")
            fh.writelines(f"{n}\t{b:.0f}\n" for n, b in skipped)


def config_datasets() -> list[str]:
    text = (ROOT / "config" / "datasets.yaml").read_text()
    return re.findall(r"^(GSE\d+):", text, flags=re.M)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("series", nargs="*", help="GSE accessions")
    p.add_argument("--all", action="store_true", help="every dataset in config/datasets.yaml")
    p.add_argument("--out", type=Path, default=ROOT / "data" / "raw")
    p.add_argument("--max-gb", type=float, default=3.0, help="skip single files larger than this")
    p.add_argument("--list-only", action="store_true")
    args = p.parse_args(argv)

    series = list(dict.fromkeys(args.series + (config_datasets() if args.all else [])))
    if not series:
        p.error("give GSE accessions or --all")
    failed = []
    for gse in series:
        try:
            fetch(gse.upper(), args.out, args.max_gb * 1e9, args.list_only)
        except Exception as exc:
            failed.append(gse)
            print(f"  [FAIL] {gse}: {exc}")
    print("\nDone." if not failed else f"\nFailed: {', '.join(failed)} — re-run to resume.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
