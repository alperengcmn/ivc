#!/usr/bin/env python3
"""Script 04 — Phase 4: depth-matched embedding and clustering.

Thin entry point; logic lives in src/ilc2tri/embed.py. Run with --help for options.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ilc2tri.embed import main

if __name__ == "__main__":
    raise SystemExit(main())
