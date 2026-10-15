#!/usr/bin/env python3
"""Script 05 — Phase 5: lineage annotation and nILC2/iILC2 subsets.

Thin entry point; logic lives in src/ilc2tri/annotate.py. Run with --help for options.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ilc2tri.annotate import main

if __name__ == "__main__":
    raise SystemExit(main())
