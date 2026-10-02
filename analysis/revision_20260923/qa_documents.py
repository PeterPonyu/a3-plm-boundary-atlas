#!/usr/bin/env python3
"""Technical PDF/render checks. Deliberately never certify visual acceptance."""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile

import fitz

HERE=Path(__file__).resolve().parent
CAP=HERE.parents[1]
MAN=CAP/"manuscript"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    documents=[]
    for stem,prefix in (("manuscript","article"),("SI","si")):
        doc=fitz.open(MAN/f"{stem}.pdf")
        beyond=[]
        blank=[]
        for i,page in enumerate(doc,1):
            if not page.get_text().strip():
                blank.append(i)
            bounds=page.rect+(-.5,-.5,.5,.5)
            for block in page.get_text("dict")["blocks"]:
                for line in block.get("lines",[]):
                    for span in line.get("spans",[]):
                        if span["text"].strip() and not bounds.contains(fitz.Rect(span["bbox"])):
                            beyond.append(dict(page=i,text=span["text"],bbox=list(span["bbox"])))
        fonts=subprocess.run(["pdffonts",str(MAN/f"{stem}.pdf")],capture_output=True,text=True,check=True).stdout
        font_lines=fonts.splitlines()[2:]
        unembedded=[line for line in font_lines if re.search(r"\sno\s+(?:yes|no)\s+(?:yes|no)\s+\d+\s+\d+\s*$",line)]
        pngs=sorted((HERE/"rendered-current").glob(prefix+"-*.png"))
        documents.append(dict(path=str((MAN/f"{stem}.pdf").relative_to(CAP)),sha256=sha(MAN/f"{stem}.pdf"),
                              pages=len(doc),rendered_pages=len(pngs),
                              all_pages_rendered=len(pngs)==len(doc),
                              render_files_nonempty=all(p.stat().st_size>10000 for p in pngs),
                              text_outside_page= beyond,blank_text_pages=blank,
                              unembedded_fonts=unembedded,type3_font="Type 3" in fonts,
                              render_sha256={p.name:sha(p) for p in pngs}))
    # The first local integration backup includes the old payload. Read it,
    # never unpack or alter it, to prove later work did not touch that bundle.
    backup=HERE/"pre-integration.08EJ4V/current-manuscript-and-contracts.tar.gz"
    historical_drift=[]
    historical_checked=0
    with tarfile.open(backup) as archive:
        for member in archive.getmembers():
            if member.isfile() and member.name.startswith("manuscript/submission/"):
                target=CAP/member.name
                before=hashlib.sha256(archive.extractfile(member).read()).hexdigest()
                historical_checked+=1
                if not target.is_file() or sha(target)!=before:
                    historical_drift.append(member.name)
    report=dict(documents=documents,render_and_bbox_pass=all(
        r["all_pages_rendered"] and r["render_files_nonempty"] and not r["text_outside_page"]
        and not r["blank_text_pages"] and not r["unembedded_fonts"] and not r["type3_font"] for r in documents),
        historical_submission_files_compared=historical_checked,historical_submission_drift=historical_drift,
        visual_review=dict(status="NOT_PERFORMED",pass_attested=False,
                           reason="functions.view_image returned: view_image is not allowed because you do not support image inputs",
                           scope="Raster rendering and bbox/font checks are not human or model visual inspection."))
    (HERE/"document-qa.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k!="documents"},indent=2))
    print(json.dumps([dict(path=r["path"],pages=r["pages"],rendered_pages=r["rendered_pages"],
                           outside_count=len(r["text_outside_page"]),sha256=r["sha256"]) for r in documents],indent=2))
    if not report["render_and_bbox_pass"] or historical_drift:
        raise SystemExit(1)


if __name__=="__main__":
    main()
