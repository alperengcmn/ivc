#!/usr/bin/env python3
"""Script 08 — Phase 8: TF (CollecTRI) and pathway (PROGENy) activity.

Thin entry point; logic lives in src/ilc2tri/tf.py. Run with --help for options.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ilc2tri.tf import main

if __name__ == "__main__":
    raise SystemExit(main())
