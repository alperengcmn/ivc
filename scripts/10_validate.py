#!/usr/bin/env python3
"""Script 10 — Phase 10: cross-dataset validation of the frozen signature (pre-registered).

Thin entry point; logic lives in src/ilc2tri/validation.py. Run with --help for options.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ilc2tri.validation import main

if __name__ == "__main__":
    raise SystemExit(main())
