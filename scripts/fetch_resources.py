#!/usr/bin/env python3
"""Download prior-knowledge resources into references/ (mouse).

    python scripts/fetch_resources.py

* MSigDB mouse collections (MH hallmark, M2 Reactome) — GMT, from the Broad release server
* CollecTRI TF→target network (mouse) — via decoupler / OmniPath
* PROGENy pathway footprints (mouse, top 500) — via decoupler / OmniPath
* UCSC mm10 refGene — gene coordinates / TSS for peak-to-gene linkage (Phase 9)
* GEO GPL16570 platform table — Affymetrix MoGene 2.0 probe → gene symbol (GSE81700)

Each file is written with a sidecar line in references/manifest.tsv (source, version, SHA-256).
"""

from __future__ import annotations

import hashlib
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "references"
MSIGDB = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/{v}/{name}.v{v}.symbols.gmt"
MSIGDB_VERSIONS = ("2025.1.Mm", "2024.1.Mm", "2023.2.Mm", "2023.1.Mm")
COLLECTIONS = {"mh.all": "msigdb_hallmark_mouse.gmt", "m2.cp.reactome": "msigdb_reactome_mouse.gmt"}


def record(path: Path, source: str, version: str) -> None:
    manifest = REF / "manifest.tsv"
    if not manifest.exists():
        manifest.write_text("file\tsource\tversion\tsha256\tdownloaded_utc\n")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with manifest.open("a") as fh:
        fh.write(f"{path.name}\t{source}\t{version}\t{digest}\t{stamp}\n")
    print(f"  [ok] {path.name}  ({version})")


def msigdb() -> bool:
    ok = True
    for name, out in COLLECTIONS.items():
        dest = REF / out
        if dest.exists():
            print(f"  [have] {out}")
            continue
        for v in MSIGDB_VERSIONS:
            url = MSIGDB.format(v=v, name=name)
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    dest.write_bytes(r.read())
                record(dest, url, v)
                break
            except Exception:
                continue
        else:
            print(f"  [FAIL] {name}: no MSigDB release reachable")
            ok = False
    return ok


MAX_BYTES = 200_000_000  # never stream more than this for a single resource
UA = {"User-Agent": "Mozilla/5.0 (ilc2-trained-immunity fetch_resources)"}
PLAIN = {
    "ucsc_mm10_refGene.txt.gz": [
        ("https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz", "mm10"),
    ],
    # Platform table only (not the multi-GB family SOFT, which embeds every sample on the platform).
    "GPL16570_table.txt": [
        (
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16570&targ=self&form=text&view=data",
            "GEO GPL16570 platform table",
        ),
    ],
}


def stream(url: str, dest: Path) -> None:
    req = urllib.request.Request(url, headers=UA)
    done = 0
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req, timeout=120) as r, tmp.open("wb") as fh:
        while chunk := r.read(1 << 20):
            done += len(chunk)
            if done > MAX_BYTES:
                raise RuntimeError(f"aborted: > {MAX_BYTES // 1_000_000} MB (wrong file?)")
            fh.write(chunk)
            print(f"\r  ... {dest.name}: {done / 1e6:.1f} MB", end="", flush=True)
    print()
    tmp.replace(dest)


def plain_files() -> bool:
    ok = True
    for name, sources in PLAIN.items():
        dest = REF / name
        if dest.exists():
            print(f"  [have] {name}")
            continue
        for url, version in sources:
            try:
                stream(url, dest)
                if name.endswith(".txt") and b"ID" not in dest.read_bytes()[:20000]:
                    dest.unlink()
                    raise RuntimeError(
                        "response is not a platform table (blocked or captcha page?)"
                    )
                record(dest, url, version)
                break
            except Exception as exc:
                print(f"  [FAIL] {name} from {url}: {type(exc).__name__}: {exc}")
        else:
            ok = False
    return ok


def _decoupler_net(kind: str):
    import decoupler as dc

    if hasattr(dc, "op"):  # decoupler >= 2.0
        fn = {"collectri": dc.op.collectri, "progeny": dc.op.progeny}[kind]
        return fn(organism="mouse"), getattr(dc, "__version__", "2.x")
    if kind == "collectri":  # decoupler 1.x
        return dc.get_collectri(organism="mouse", split_complexes=False), dc.__version__
    return dc.get_progeny(organism="mouse", top=500), dc.__version__


def networks() -> bool:
    ok = True
    for kind in ("collectri", "progeny"):
        dest = REF / f"{kind}_mouse.tsv"
        if dest.exists():
            print(f"  [have] {dest.name}")
            continue
        try:
            net, version = _decoupler_net(kind)
            net.to_csv(dest, sep="\t", index=False)
            record(dest, f"decoupler/OmniPath {kind}", f"decoupler {version}")
        except Exception as exc:
            print(f"  [FAIL] {kind}: {type(exc).__name__}: {exc}")
            ok = False
    return ok


def main() -> int:
    REF.mkdir(exist_ok=True)
    print("== MSigDB")
    a = msigdb()
    print("== TF / pathway networks")
    b = networks()
    print("== Genome annotation / array platforms")
    c = plain_files()
    return 0 if a and b and c else 1


if __name__ == "__main__":
    sys.exit(main())
