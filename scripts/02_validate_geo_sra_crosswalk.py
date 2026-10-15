#!/usr/bin/env python3
"""Script 02 — GEO → BioSample → SRA Experiment → SRR crosswalk validation.

Thin entry point so the script runs from a fresh clone without installing the
package. Logic lives in src/ilc2tri/validate_crosswalk.py. Run with --help for options.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ilc2tri.validate_crosswalk import main

if __name__ == "__main__":
    raise SystemExit(main())
