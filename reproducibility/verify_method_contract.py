"""Check bound source identity; this does not attest scientific validity."""
from pathlib import Path
import hashlib,json
CAP=Path(__file__).resolve().parents[1]
record=json.loads((CAP/'reproducibility/method-contract.json').read_text())
fail=[]
for item in record['locked_files']:
    p=CAP/item['path']
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=item['sha256']:
        fail.append(item['path'])
if fail:raise SystemExit('SOURCE IDENTITY CHANGED; review protocol before rebuilding: '+', '.join(fail))
print(f"{record['paper']}: {len(record['locked_files'])} bound inputs/code files match. Status: {record['status']}")
