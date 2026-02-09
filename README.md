# Trove Raman Fractions Analysis

Analysis code and data for the Trove Raman Fractions paper.

## Overview

This repository contains the Jupyter notebook analysis for comparing SERS (Surface-Enhanced Raman Spectroscopy) and Raman spectra from Trove fraction samples.

## Repository Structure

```
trove-raman-paper/
├── notebooks/
│   └── trove-raman-fractions.ipynb  # Main analysis notebook
├── raman_analysis/                   # Python package for Raman data analysis
│   ├── __init__.py
│   └── utils.py                      # Utilities for parsing and plotting spectra
├── data/                             # Spectroscopy data files
│   ├── 2025-03-18-ticks/            # Raman data
│   └── Trove-SERS-attempt3/         # SERS data
│       └── 2025-07-24/
├── pyproject.toml                    # Package dependencies
└── README.md
```

## Setup

### Prerequisites

- Python >= 3.12
- pip or uv package manager

### Installation

1. Clone this repository:
```bash
git clone <repository-url>
cd trove-raman-paper
```

2. Install dependencies:
```bash
pip install -e .
```

Or using uv:
```bash
uv pip install -e .
```

### Dependencies

The main dependencies are:
- `pandas` - Data manipulation
- `plotly` - Interactive plotting
- `ramanspy` - Raman spectroscopy preprocessing
- `numpy` - Numerical computing
- `arcadia_pycolor` - Arcadia color schemes for plotting

## Usage

Open and run the Jupyter notebook:

```bash
jupyter lab notebooks/trove-raman-fractions.ipynb
```

The notebook performs:
1. Loading Raman and SERS data from multiple experiments
2. Preprocessing using RamanSPy pipeline (cropping, despiking, denoising, baseline correction)
3. Visualization and comparison of spectra

## Data

The data directory contains:
- **2025-03-18-ticks/**: Raman spectra from whole extract measurements
- **Trove-SERS-attempt3/2025-07-24/**: SERS spectra from fraction measurements

Each data directory includes:
- Spectroscopy data files (.txt, .csv, or .tsv)
- `sample_table.csv` with metadata mapping file IDs to sample names and groups

## Citation

If you use this code or data, please cite:

[Add citation information here]

## License

MIT License
