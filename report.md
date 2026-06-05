# Assignment 1: Characterising a SAMBA AMOC Time Series

BBF1004 Yuxiang Luan

## Data choice

This analysis uses the SAMBA 34.5S observing array from AMOCatlas and the
`UPPER_TRANSPORT` variable. This series is the upper-cell volume transport
anomaly in Sverdrups (Sv). I chose SAMBA because it provides an AMOC-related
time series from the South Atlantic and does not duplicate the MOCHA/RAPID 26N
dataset assigned to another student.

The record spans 2013-09-12 to 2017-07-16. It contains 1,404 daily samples, with
a median sampling interval of 1.0 day. There are no missing values in the
standardised series, so no interpolation was needed. The gap-filling step in the
script is kept for reproducibility, but it leaves this particular series
unchanged.

## Time-domain statistics

The upper-cell transport anomaly has a mean of -0.003 Sv, which is essentially
zero as expected for an anomaly series. The standard deviation is 15.51 Sv,
showing substantial variability around the mean. The median is 0.35 Sv. The
minimum and maximum are -50.28 Sv and 52.69 Sv, giving a total range of
102.97 Sv. The histogram shows that most values are clustered near zero, but the
series also includes large positive and negative transport anomalies.

![SAMBA upper-cell transport distribution](figures/assignment1_distribution_histogram.png)

The raw time series contains strong day-to-day and multi-week variability. I
applied a 90-day Tukey-window low-pass filter using a centred 91-sample rolling
window. The filtered curve removes much of the short-timescale variability and
highlights slower changes in the upper-cell transport anomaly.

![Raw and filtered SAMBA upper-cell transport](figures/assignment1_timeseries_raw_filtered.png)

## Frequency-domain statistics

I estimated the power spectral density using Welch's method with Hann windows,
1024-sample segments, and 50 percent overlap. Since the data are daily, each
segment spans 1024 days. This segment length is long enough to retain useful
low-frequency resolution for the relatively short SAMBA record, while still
using overlapped segment averaging rather than a single raw periodogram.

The integrated Welch spectrum satisfies the Parseval variance check. The ratio
between integrated PSD and time-domain variance is 1.002 for the raw series,
which is very close to 1. This indicates that the spectral normalisation is
consistent with the variance of the time series. The largest spectral peak in
this estimate corresponds to a timescale of about 256 days. Overall, the
spectrum has stronger low-frequency than high-frequency variance, suggesting a
red spectrum with variability concentrated at multi-month and longer
timescales.

The 90-day Tukey low-pass filter reduces spectral power at frequencies above
about 1/90 cycles per day. The filtered spectrum therefore lies below the raw
spectrum at higher frequencies, while retaining more of the lower-frequency
variance. This behaviour is consistent with the goal of isolating slower AMOC
transport variability.

![Raw and filtered SAMBA upper-cell spectrum](figures/assignment1_spectrum_raw_filtered.png)

## Reproducibility

The analysis can be reproduced with:

```bash
python assignment_analysis.py
pytest -q
```

The script downloads/loads the SAMBA dataset through AMOCatlas, writes the
figures to `figures/`, and writes the numerical summary to
`outputs/assignment1_summary.json`.
