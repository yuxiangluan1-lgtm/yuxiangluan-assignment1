# Assignment 1: Characterising a RAPID AMOC Component Time Series

## Data choice

This analysis uses the RAPID 26N AMOC observing array and the `t_ek10` variable
from `moc_transports.nc`. This variable is the Ekman transport component, in
Sverdrups (Sv). I chose it because it is a physically meaningful AMOC component
but is not the prohibited headline overturning series `moc_mar_hc10`.

The record spans 2004-04-02 to 2023-02-11. The median sampling interval is
0.5 days, so the regular grid is effectively a 12-hour time series. The series
contains 13,779 samples, including 20 missing values. The longest missing run is
10 samples, or 5.0 days. For the spectral calculations I filled these gaps using
linear interpolation, keeping the original gappy series for the descriptive
statistics.

## Time-domain statistics

The Ekman transport has a mean of 3.77 Sv, a standard deviation of 3.44 Sv, and a
median of 3.88 Sv. Its minimum and maximum are -13.00 Sv and 18.29 Sv,
respectively, giving a range of 31.29 Sv. The distribution is therefore centred
on a small positive transport but includes substantial positive and negative
events. Compared with the mean, the standard deviation is large, which indicates
that high-frequency wind-driven variability is a major part of this component.

The time series plot shows rapid synoptic-to-seasonal variability throughout the
record. I applied a 90-day Tukey-window low-pass filter with a centred rolling
window of 181 samples. The filtered curve suppresses short-period variability and
shows slower changes more clearly, while preserving the broad multi-month to
interannual structure.

![Raw and filtered Ekman transport](figures/assignment1_timeseries_raw_filtered.png)

## Frequency-domain statistics

I estimated the power spectral density using Welch's method with Hann windows,
1024-sample segments, and 50 percent overlap. With a 0.5-day sampling interval,
each segment spans 512 days. This choice averages several independent segments,
reducing the noise of a raw periodogram while retaining enough resolution to
separate synoptic, seasonal, and longer timescales.

The integrated Welch spectrum satisfies the Parseval variance check: the ratio
between integrated PSD and time-domain variance is 0.996 for the raw filled
series. This is close to 1, so the PSD normalization is consistent with the
variance of the time series. The largest Welch spectral peak occurs at an
approximately 512-day timescale. The spectrum has more energy at low frequencies
than at the highest resolved frequencies, which is consistent with a redder
low-frequency component superimposed on strong higher-frequency wind-driven
variability.

The 90-day Tukey low-pass filter reduces spectral power at frequencies above
about 1/90 cycles per day. In the filtered spectrum, the high-frequency part is
strongly attenuated, while the lower-frequency variance is retained. This is the
expected behaviour of a low-pass filter and makes the slower AMOC-component
variability easier to interpret.

![Raw and filtered Ekman spectrum](figures/assignment1_spectrum_raw_filtered.png)

## Reproducibility

The analysis can be reproduced with:

```bash
python assignment_analysis.py
pytest -q
```

The script writes the figures to `figures/` and the numerical summary to
`outputs/assignment1_summary.json`.
