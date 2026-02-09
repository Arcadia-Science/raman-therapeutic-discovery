"""
Raman Analysis Package

A Python package for analyzing Raman spectroscopy data with support for multiple instruments
including Horiba, Wasatch, and OpenRaman spectrometers.
"""

from .utils import (
    collect_spectra_from_dir,
    parse_horiba_output,
    parse_openraman_output,
    parse_wasatch_output,
    plot_spectra_in_dir,
)

__version__ = "0.1.0"
__all__ = [
    "collect_spectra_from_dir",
    "parse_horiba_output", 
    "parse_openraman_output",
    "parse_wasatch_output",
    "plot_spectra_in_dir",
    "sort_files_to_plot",
]