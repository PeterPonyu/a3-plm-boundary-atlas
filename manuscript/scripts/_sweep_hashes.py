#!/usr/bin/env python3
"""Throwaway audit: compare every recorded hash/page field to measured bytes."""
import hashlib
import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
FLAT = ROOT / "submission/flat_Proteins"


def sha(p):
    p = pathlib.Path(p)
    if not p.exists():
        return None
    return hashlib.sha256(p.read_bytes()).hexdigest()


def pages(p):
    p = pathlib.Path(p)
    if not p.exists():
        return None
    out = subprocess.run(["pdfinfo", str(p)], capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith("Pages:"):
            return int(line.split()[1])
    return None


findings = []


def check(label, recorded, measured):
    status = "OK " if recorded == measured else "MISMATCH"
    if recorded != measured:
        findings.append((label, recorded, measured))
    print(f"{status}  {label}\n          recorded={recorded}\n          measured={measured}")


spec = json.loads((ROOT / "release-spec.json").read_text())
att = json.loads((ROOT / "attestation/build-attestation.json").read_text())
mirror = json.loads((FLAT / "release-spec.json").read_text())

print("=" * 78)
print("SOURCE release-spec.json")
print("=" * 78)
main = spec["documents"]["main"]
check("documents.main.sha256 -> " + main["pdf"], main["sha256"], sha(ROOT / main["pdf"]))
check("documents.main.pages", main["pages"], pages(ROOT / main["pdf"]))
if main.get("pdf_sha256") is not None:
    check("documents.main.pdf_sha256", main["pdf_sha256"], sha(ROOT / main["pdf"]))
else:
    print("NOTE      documents.main.pdf_sha256 is null (unbound, not stale)")

for e in spec["package"]["upload_payload"]:
    check("package.upload_payload " + e["path"], e["sha256"], sha(FLAT / e["path"]))

for e in spec.get("standalone_assets") or []:
    check("standalone_assets " + e["path"], e["sha256"], sha(ROOT / e["path"]))

drift = 0
for e in spec["build_inputs"]:
    m = sha(ROOT / e["path"])
    if m != e["sha256"]:
        drift += 1
        findings.append(("build_inputs " + e["path"], e["sha256"], m))
        print(f"MISMATCH  build_inputs {e['path']}\n          recorded={e['sha256']}\n          measured={m}")
print(f"NOTE      build_inputs: {len(spec['build_inputs'])} entries, {drift} drifted")

print()
print("=" * 78)
print("attestation/build-attestation.json")
print("=" * 78)
check("builder.sha256", att["builder"]["sha256"], sha(ROOT / "scripts" / att["builder"]["path"]))
check("source_release_spec.sha256", att["source_release_spec"]["sha256"], sha(ROOT / "release-spec.json"))
check("package_release_spec.sha256", att["package_release_spec"]["sha256"], sha(FLAT / "release-spec.json"))
check("document.sha256", att["document"]["sha256"], sha(FLAT / "manuscript.pdf"))
check("document.pages", att["document"]["pages"], pages(FLAT / "manuscript.pdf"))
check("log.sha256", att["log"]["sha256"], sha(ROOT / att["log"]["path"]))
check("log_policy.sha256", att["log_policy"]["sha256"], sha(ROOT / att["log_policy"]["path"]))
for k, v in att["upload_payload"].items():
    check("att.upload_payload " + k, v, sha(FLAT / k))
check("att.build_inputs_digest", att["build_inputs_digest"], spec["build_inputs_digest"])

print()
print("=" * 78)
print("payload mirrored release-spec.json")
print("=" * 78)
for d in mirror["documents"]:
    check("mirror documents[].sha256", d["sha256"], sha(FLAT / d["pdf"]))
    check("mirror documents[].pages", d["pages"], pages(FLAT / d["pdf"]))
for k, v in (mirror.get("standalone_assets") or {}).items():
    check("mirror standalone_assets " + k, v, sha(FLAT / k))
for e in mirror["upload_payload"]:
    check("mirror upload_payload " + e["path"], e["sha256"], sha(FLAT / e["path"]))
check("mirror build_inputs_digest", mirror["build_inputs_digest"], spec["build_inputs_digest"])

print()
print("=" * 78)
print("SHA256SUMS.txt")
print("=" * 78)
listed = set()
for line in (FLAT / "SHA256SUMS.txt").read_text().splitlines():
    if not line.strip():
        continue
    h, name = line.split(None, 1)
    name = name.strip()
    listed.add(name)
    check("SHA256SUMS " + name, h, sha(FLAT / name))
present = {p.name for p in FLAT.iterdir() if p.is_file()}
print("NOTE      files in payload dir NOT in SHA256SUMS:", sorted(present - listed - {"SHA256SUMS.txt"}))
print("NOTE      files in SHA256SUMS not on disk:", sorted(listed - present))

print()
print("=" * 78)
print("visual-registry.json (staged + payload)")
print("=" * 78)
for cand in [FLAT / "visual-registry.json", ROOT / "visual-review-source.json"]:
    if not cand.exists():
        continue
    vr = json.loads(cand.read_text())
    print(f"--- {cand.relative_to(ROOT)} top keys: {list(vr.keys()) if isinstance(vr, dict) else type(vr)}")
    blob = json.dumps(vr)
    import re
    for m in re.finditer(r'"(sha256|pages|page|page_count)"\s*:\s*("?[0-9a-f]{16,}"?|\d+)', blob):
        print("    field:", m.group(0)[:120])

print()
print("=" * 78)
print(f"TOTAL MISMATCHES: {len(findings)}")
for label, r, m in findings:
    print(f"  - {label}: recorded {r} != measured {m}")
sys.exit(0)
