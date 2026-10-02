#!/usr/bin/env python3
"""Check the current A3 source/PDF boundary and record layout diagnostics."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

import fitz

HERE = Path(__file__).resolve().parent
CAP = HERE.parents[1]
MAN = CAP/"manuscript"


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    source = (MAN/"manuscript.tex").read_text()
    clinical = (MAN/"text/clinical_revision.tex").read_text()
    methods = (MAN/"text/methods.tex").read_text()
    text = "\n".join(p.get_text() for p in fitz.open(MAN/"manuscript.pdf"))
    tests = {
        "original_source_includes_revised_clinical_sections": r"\input{text/clinical_revision}" in source,
        "no_old_unified_calibration_headline": all(x not in text for x in ("7/20 pattern", "0.2487", "0.0204", "0.0131")),
        "no_unmeasured_transfer_failure": "does not transfer to ClinVar" not in source,
        "historical_protocol_distinguished": "seed 1" in methods and "all 1,093" in methods,
        "pdf_same_score_clinical_recovery": "4,677" in text and "0.7841" in text and "0.1412" in text,
        "new_intervals_have_explicit_generating_protocol": "2,000" in methods and "1,000 refitting" in methods and "multiplicity" in methods,
        "no_false_still_missing_scores": "row-level recovery remains incomplete" not in text and "score-version audit remains open" not in text,
        "new_ties_and_identity_scope_disclosed": "equal scores" in methods and "does not establish independence" in methods,
        "hold_on_pdf_first_page": "not submission-ready" in fitz.open(MAN/"manuscript.pdf")[0].get_text(),
        "hold_on_si_first_page": "not submission-ready" in fitz.open(MAN/"SI.pdf")[0].get_text(),
        "method_contract_checks_pass": subprocess.run(["python",str(CAP/"reproducibility/verify_method_contract.py")],capture_output=True).returncode==0,
        "clinical_numerical_checks_pass": subprocess.run(["python",str(HERE/"validate_clinical.py")],capture_output=True).returncode==0,
    }
    for stem in ("manuscript", "SI"):
        log=(MAN/f"{stem}.log").read_text(errors="replace")
        tests[stem+"_no_overfull_or_undefined"] = not re.search(r"Overfull|undefined|LaTeX Warning",log)
    release=json.loads((MAN/"release-current.json").read_text())
    drift=[r["path"] for r in release["build_inputs"]+release["documents"] if digest(MAN/r["path"])!=r["sha256"]]
    tests["release_binds_current_sources_and_pdfs"] = not drift
    tests["no_submission_ready_attestation"] = release["scientific_ready_attested"] is False
    body="\n".join((MAN/p).read_text() for p in ("manuscript.tex","text/clinical_revision.tex","text/methods.tex","text/abstract.tex"))
    tests["generated_clinical_macros_connected"] = r"\input{text/clinical_values_20260923}" in body
    record=dict(checks=tests,all_pass=all(tests.values()),drift=drift,
                visual_review="Machine-only verification; separately hash-bound rendered-page inspection is required.")
    (HERE/"integration-validation.json").write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps(record,indent=2))
    if not record["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
