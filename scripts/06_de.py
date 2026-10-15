#!/usr/bin/env python3
"""Script 06 — Phase 6: exploratory library-consistency differential expression in ILC2s.

Thin entry point; logic lives in src/ilc2tri/de.py. Run with --help for options.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ilc2tri.de import main

if __name__ == "__main__":
    raise SystemExit(main())
