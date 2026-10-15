#!/usr/bin/env python3
"""Script 01 — structural validation of GSE295590 sample metadata.

Thin entry point so the script runs from a fresh clone without installing the
package. Logic lives in src/ilc2tri/validate_samples.py. Run with --help for options.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ilc2tri.validate_samples import main

if __name__ == "__main__":
    raise SystemExit(main())
