#!/usr/bin/env python3
"""Make Cairo's fixed-width PDF creation-date metadata reproducible."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import fitz


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: normalize_a3_pdf_metadata.py PDF")
    pdf = Path(sys.argv[1])
    document = fitz.open(pdf)
    metadata = document.metadata
    fixed_date = "D:19700101000000+00'00'"
    metadata["creationDate"] = fixed_date
    metadata["modDate"] = fixed_date
    document.set_metadata(metadata)
    normalized = pdf.with_suffix(".normalized.pdf")
    document.save(normalized, garbage=4, deflate=True, no_new_id=True)
    document.close()
    os.replace(normalized, pdf)


if __name__ == "__main__":
    main()
