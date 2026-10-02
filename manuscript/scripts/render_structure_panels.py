#!/usr/bin/env python3
"""Optional coordinate rerender; the default build uses bundled verified PNGs."""
import os,shutil,subprocess
from pathlib import Path
root=Path(__file__).resolve().parents[1]
exe=os.environ.get('PYMOL_EXE') or shutil.which('pymol')
if not exe:raise SystemExit('Set PYMOL_EXE to an existing PyMOL executable.')
for p in sorted((root/'objects').rglob('render.pml')):
    subprocess.run([exe,'-cq',p.name],cwd=p.parent,check=True)
