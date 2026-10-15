#!/usr/bin/env python3
"""Script 09 — Phase 9: GSE172258 ATAC-seq differential accessibility vs the frozen signature.

Thin entry point; logic lives in src/ilc2tri/chromatin.py. Run with --help for options.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ilc2tri.chromatin import main

if __name__ == "__main__":
    raise SystemExit(main())
