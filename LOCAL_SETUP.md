# Local installation on Max's Mac

This repository bundle was generated in ChatGPT's sandbox and cannot be written directly into the host Mac filesystem. After downloading/extracting the bundle, place the folder at:

```text
/Users/maxpiper/CRD_DAP
```

Then, from Terminal:

```bash
cd /Users/maxpiper/CRD_DAP

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -e ".[science,dev]"

python -m pytest -q
```

The current infrastructure tests should pass before Script 1 development begins.

When ready to initialize Git locally:

```bash
git init
git add .
git commit -m "Initialize CRD_DAP pipeline architecture"
```

Do not commit large FITS cubes, master arcs, XSL libraries, PyMorph VAC files, or generated run products. The included `.gitignore` is already set up to exclude the normal large-data/product locations and common array/FITS outputs.
