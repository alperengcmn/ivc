#!/usr/bin/env python3
"""Script 03 — Phase 3 single-cell QC for GSE295590 (report-only until thresholds are agreed).

Thin entry point; logic lives in src/ilc2tri/qc.py. Run with --help for options.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ilc2tri.qc import main

if __name__ == "__main__":
    raise SystemExit(main())
