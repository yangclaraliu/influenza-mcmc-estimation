# Influenza ABC-MCMC estimation

This repository estimates household and school influenza transmission
parameters for elementary, middle, and high schools using ABC-MCMC.
Kindergarten is intentionally excluded.

## Requirements

- Python 3.10 or newer
- Authorized local copies of the four input workbooks described in
  [DATA_FORMAT.md](DATA_FORMAT.md)

Create an environment and install the dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Place the input workbooks in the repository root. They are intentionally
excluded from Git because they contain individual-level health information.

## Run

Perform a short smoke test:

```bash
python "MCMC estimation(all schools).py" --quick
```

Run the full analysis with the default settings:

```bash
python "MCMC estimation(all schools).py"
```

Estimate selected school types:

```bash
python "MCMC estimation(all schools).py" --schools elementary middle
```

Use `--help` to see all sampling and calibration options. Output CSV and Excel
files are written to the repository root and ignored by Git.

## Data availability

The original input data are not publicly distributed because they include
individual-level health and demographic fields. Researchers with appropriate
authorization can prepare compatible local files using
[DATA_FORMAT.md](DATA_FORMAT.md).
