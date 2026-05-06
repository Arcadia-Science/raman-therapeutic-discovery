# Raman Spectroscopy and SERS of Tick Lipids

This repo can be used to generate the figure in the pub [**Raman spectroscopy may assist natural therapeutic discovery**](https://doi.org/10.57844/arcadia-2mmy-8h28). Open and run `notebooks/generate_figure.ipynb`. It also contains the raw data from the experiments.


## Setup

### Prerequisites

- [uv](https://docs.astral.sh/uv/) package manager
- Python 3.12

### Usage

1. Clone this repository:
```bash
git clone https://github.com/Arcadia-Science/raman-therapeutic-discovery/
cd raman-bioprospecting
```

2. Install dependencies using uv:
```bash
uv sync
```

This will create a virtual environment and install all dependencies specified in `pyproject.toml`.

Navigate to `notebooks/` and run the Jupyter notebook `generate_figure.ipynb`.

## Raw data

The `data` directory contains:
- **whole-tick-extract/**: Raman spectra from whole tick extract measurements (Round 1)
- **sers-on-fraction/**: SERS spectra from tick fraction measurements (Round 2)

Each data subdirectory includes:
- Raw spectroscopy data files (.csv) containing acquisition metadata and spectra
- `sample_table.csv` metadata table
