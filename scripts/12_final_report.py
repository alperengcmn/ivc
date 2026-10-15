#!/usr/bin/env python3
"""Script 12 — Phase 12: assemble reports/FINAL_REPORT.md from result files.

Thin entry point; logic lives in src/ilc2tri/final_report.py.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ilc2tri.final_report import main

if __name__ == "__main__":
    raise SystemExit(main())
