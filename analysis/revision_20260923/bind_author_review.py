#!/usr/bin/env python3
"""Bind the actual A3 author-review source/PDF pair, never certify readiness.

Run only after successful figure/text checks and compilation. Earlier release
and review records remain in the pre-integration backup and historical payload.
"""
from pathlib import Path
import argparse
import hashlib
import json

import fitz

HERE = Path(__file__).resolve().parent
CAP = HERE.parents[1]
MAN = CAP / "manuscript"
STATUS = "AUTHOR_REVIEW_SAME_SCORE_INTEGRATED_NOT_SUBMISSION_READY"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lock-inputs", action="store_true",
                    help="Accept intentional source changes only; do not attest existing PDFs")
    args = ap.parse_args()
    required = ["manuscript.tex", "SI.tex", "text/clinical_revision.tex", "text/methods.tex", "text/abstract.tex"]
    for name in required:
        if not (MAN/name).is_file():
            raise SystemExit(f"Required manuscript source missing: {name}")
    source = (MAN/"text/clinical_revision.tex").read_text()
    if "row-level" not in source or "historical" not in source.lower():
        raise SystemExit("Clinical source-version disclosure missing")
    summary = json.loads((HERE/"clinical-unified/summary.json").read_text())
    if digest(HERE/"row-scores/scores.csv") != summary["input_sha256"]:
        raise SystemExit("Clinical score/analysis binding failed")
    if digest(HERE/"analyze_clinical.py") != summary["analysis_sha256"]:
        raise SystemExit("Clinical analysis was not produced by current code")
    inputs = set()
    for pattern in ("*.tex", "text/*.tex", "captions/*.tex", "figures/*.pdf", "figures/*.png",
                    "figures/plotdata/*.csv", "figures/plotdata/*.json", "scripts/R/*.R",
                    "scripts/*.py", "scripts/*.sh", "scripts/tests/*.py"):
        inputs.update(p for p in MAN.glob(pattern) if p.is_file())
    bound = [CAP / item["path"] for item in json.loads((CAP/"reproducibility/method-contract.json").read_text())["locked_files"]]
    bound += [p for p in inputs if p.suffix not in (".pdf", ".png")]
    for pattern in ("*.py", "row-scores/*.csv", "row-scores/*.json", "clinical-unified/*"):
        bound += [p for p in HERE.glob(pattern) if p.is_file()]
    bound += [CAP/"REPRODUCIBILITY-CURRENT.md", CAP/"reproducibility/METHOD-AUDIT.md",
              CAP/"analysis/revision_20260922/candidate-v1-checkpoint/download-provenance.json",
              MAN/"FIGURE-PROVENANCE.md"]
    contract = dict(paper="A3",status=STATUS,scientific_ready_attested=False,
                    scope=["4677 newly inferred variant scores pass 24 N6c full-sequence/window gates; all current clinical calculations use that score version.",
                           "Per-gene AUROC intervals newly generated; no old calibration intervals, degradation counts or clinical deployment conclusions promoted.",
                           "Tie-aggregated isotonic and identity-weighted refits remain within-panel and conditional on gene/class composition.",
                           "gnomAD and historical predictor sensitivity are summary-only. Source/PDF hashes attest identity, not readiness."],
                    locked_files=[dict(path=str(p.relative_to(CAP)),sha256=digest(p)) for p in sorted(set(bound))])
    if args.lock_inputs:
        (CAP/"reproducibility/method-contract.json").write_text(json.dumps(contract,indent=2)+"\n")
        print(f"Accepted {len(set(bound))} current input hashes; PDFs have NOT been rebound")
        return
    previous = json.loads((CAP/"reproducibility/method-contract.json").read_text())
    drift=[item["path"] for item in previous["locked_files"] if digest(CAP/item["path"])!=item["sha256"]]
    if drift:
        raise SystemExit("Unaccepted input drift; review before --lock-inputs: "+", ".join(drift))
    documents = []
    for stem in ("manuscript", "SI"):
        pdf = MAN/f"{stem}.pdf"
        doc = fitz.open(pdf)
        text = "\n".join(p.get_text() for p in doc)
        if "not submission-ready" not in text:
            raise SystemExit(f"Unmarked or stale PDF: {pdf}")
        documents.append(dict(path=pdf.name, pages=len(doc), sha256=digest(pdf)))
    record = dict(schema_version="A3-author-review-20260923", entrypoint="scripts/rebuild_current.sh",
                  status=STATUS, scientific_ready_attested=False,
                  scope="Original near-ready manuscript revised in place. Nonclinical evidence retained; same-score clinical discrimination, OOF calibration and identity-grouped refitting integrated. Historical gnomAD sensitivity remains separately labelled.",
                  clinical_input_sha256=summary["input_sha256"],
                  clinical_protocol=summary["protocol"],
                  unresolved=["Separate visual acceptance record required; machine build cannot attest it", "Independent scientific review of new analysis", "No external clinical validation", "Antibody interval source identities unresolved; archived summaries kept separate without significance inference"],
                  build_inputs=[dict(path=str(p.relative_to(MAN)), sha256=digest(p)) for p in sorted(inputs)],
                  documents=documents,
                  historical_payload="submission/flat_Proteins: historical only; not this revision and not upload authority")
    for name in ("release-current.json", "release-spec.json"):
        (MAN/name).write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(dict(status=STATUS,documents=documents,manuscript_bound_inputs=len(inputs),method_bound_files=len(set(bound))),indent=2))


if __name__ == "__main__":
    main()
