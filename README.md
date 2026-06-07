# Assignment 1: NAC AMOC Time-Series Characterisation

This repository contains my solution for Assignment 1 in Data Analysis in
Physical Oceanography. The analysis characterises the North Atlantic Current
transport proxy (`TRANS_NAC_PROXY`) in the time and frequency domains.

## Contents

- `assignment_analysis.py` runs the full analysis and writes figures plus a JSON
  numerical summary.
- `spectra_filtering/` contains the time-domain, spectral, filtering, and data
  helper functions used by the analysis.
- `tests/test_assignment1_requirements.py` contains the required tests for the
  mean calculation and Parseval PSD check.
- `figures/assignment1_timeseries_raw_filtered.png` shows the raw and 5-year
  Tukey-filtered time series.
- `figures/assignment1_spectrum_raw_filtered.png` shows the raw and filtered
  Welch spectra.
- `figures/assignment1_distribution_histogram.png` shows the distribution of
  the selected NAC transport proxy series.
- `figures/assignment1_filter_response.png` compares the 5-year Tukey filter
  response with a same-width boxcar filter response.
- `report.md` is the short written report.

## Reproduce

```bash
pip install -r requirements.txt
python assignment_analysis.py
pytest -q
```

The data are loaded with AMOCatlas using `read.nac()`.
