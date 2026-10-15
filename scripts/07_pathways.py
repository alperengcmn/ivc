#!/usr/bin/env python3
"""Script 07 — Phase 7: pathway enrichment and frozen discovery signature.

Thin entry point; logic lives in src/ilc2tri/pathways.py. Run with --help for options.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ilc2tri.pathways import main

if __name__ == "__main__":
    raise SystemExit(main())
