# Raman Spectroscopy and SERS of Tick Lipids

This repo can be used to generate the figure in the pub **Raman spectroscopy-assisted natural therapeutic discovery**. Open and run `notebooks/generate_figure.ipynb`


## Setup

### Prerequisites

- [uv](https://docs.astral.sh/uv/) package manager
- Python 3.12.11

### Usage

1. Clone this repository:
```bash
git clone <repository-url>
cd 2026-rs-tick-lipids
```

2. Install dependencies using uv:
```bash
uv sync
```

This will create a virtual environment and install all dependencies specified in `pyproject.toml`.

3. Activate the virtual environment and install the `raman_analysis` package in editable mode:
```bash
source .venv/bin/activate  # On macOS/Linux
# or
.venv\Scripts\activate  # On Windows

pip install -e .
```

Open and run the Jupyter notebook `notebooks/generate_figure.ipynb`

## Raw data

The data directory contains:
- **2025-03-18-ticks/**: Raman spectra from tick extract measurements
- **Trove-SERS-attempt3/2025-07-24/**: SERS spectra from tick fraction measurements

Each data directory includes:
- Spectroscopy data files (.txt, .csv, or .tsv)
- `sample_table.csv` with metadata mapping file IDs to sample names and groups
