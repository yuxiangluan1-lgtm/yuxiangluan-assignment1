"""Run Assignment 1 analysis for a SAMBA AMOC time series."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from amocatlas import read

from spectra_filtering.analysis import summary_stats
from spectra_filtering.data_io import fill_gaps
from spectra_filtering.filters import tukey_lowpass
from spectra_filtering.spectra import parseval_ratio, welch_psd


FIGURE_DIR = Path("figures")
OUTPUT_DIR = Path("outputs")

ARRAY_NAME = "SAMBA 34.5S"
VARIABLE = "UPPER_TRANSPORT"
VARIABLE_LABEL = "SAMBA 34.5S upper-cell transport anomaly"
UNITS = "Sv"
SEGMENT_LENGTH = 1024
OVERLAP = 0.5
LOWPASS_DAYS = 90.0


def _longest_nan_run(values: np.ndarray) -> int:
    missing = ~np.isfinite(values)
    longest = run = 0
    for item in missing:
        if item:
            run += 1
            longest = max(longest, run)
        else:
            run = 0
    return longest


def _dominant_timescale_days(freq: np.ndarray, psd: np.ndarray) -> float:
    valid = freq > 0
    peak = int(np.argmax(psd[valid]))
    return float(1.0 / freq[valid][peak])


def _load_series() -> tuple[np.ndarray, np.ndarray, float]:
    ds = read.samba()
    da = ds[VARIABLE]
    time = da["TIME"].values
    values = da.values.astype("float64")
    dt_days = float(np.median(np.diff(time)) / np.timedelta64(1, "D"))
    return time, values, dt_days


def run() -> dict[str, object]:
    FIGURE_DIR.mkdir(exist_ok=True)
    OUTPUT_DIR.mkdir(exist_ok=True)

    time, values, dt_days = _load_series()
    filled = fill_gaps(values)
    lp_window = int(round(LOWPASS_DAYS / dt_days))
    if lp_window % 2 == 0:
        lp_window += 1
    filtered = tukey_lowpass(filled, window=lp_window)

    freq, psd = welch_psd(
        filled,
        dt_days=dt_days,
        segment_length=SEGMENT_LENGTH,
        overlap=OVERLAP,
    )
    freq_f, psd_f = welch_psd(
        filtered,
        dt_days=dt_days,
        segment_length=SEGMENT_LENGTH,
        overlap=OVERLAP,
    )

    stats = summary_stats(values)
    ratio = parseval_ratio(filled, freq, psd)
    ratio_filtered = parseval_ratio(filtered, freq_f, psd_f)

    time_start = np.datetime_as_string(time[0], unit="D")
    time_end = np.datetime_as_string(time[-1], unit="D")
    summary = {
        "array": ARRAY_NAME,
        "variable": VARIABLE,
        "variable_label": VARIABLE_LABEL,
        "time_start": time_start,
        "time_end": time_end,
        "dt_days": dt_days,
        "lowpass_days": LOWPASS_DAYS,
        "lowpass_window_samples": lp_window,
        "welch_segment_length": SEGMENT_LENGTH,
        "welch_overlap": OVERLAP,
        "longest_gap_samples": _longest_nan_run(values),
        "longest_gap_days": _longest_nan_run(values) * dt_days,
        "dominant_timescale_days": _dominant_timescale_days(freq, psd),
        "parseval_ratio": ratio,
        "parseval_ratio_filtered": ratio_filtered,
        "stats": stats,
    }

    with (OUTPUT_DIR / "assignment1_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    fig, ax = plt.subplots(figsize=(10, 4.8), constrained_layout=True)
    ax.plot(time, filled, color="#31688e", linewidth=0.8, alpha=0.7, label="Gap-filled raw")
    ax.plot(time, filtered, color="#b7372f", linewidth=1.6, label=f"{LOWPASS_DAYS:.0f}-day Tukey low-pass")
    ax.set_title(f"{VARIABLE_LABEL} time series")
    ax.set_ylabel(f"Transport ({UNITS})")
    ax.set_xlabel("Time")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(FIGURE_DIR / "assignment1_timeseries_raw_filtered.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    ax.hist(values[np.isfinite(values)], bins=32, color="#31688e", alpha=0.82, edgecolor="white")
    ax.axvline(stats["mean"], color="#b7372f", linewidth=1.6, label=f"Mean = {stats['mean']:.2f} {UNITS}")
    ax.set_title(f"{VARIABLE_LABEL} distribution")
    ax.set_xlabel(f"Transport anomaly ({UNITS})")
    ax.set_ylabel("Count")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(FIGURE_DIR / "assignment1_distribution_histogram.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    positive = freq > 0
    positive_f = freq_f > 0
    ax.loglog(freq[positive], psd[positive], color="#31688e", linewidth=1.1, label="Raw Welch PSD")
    ax.loglog(freq_f[positive_f], psd_f[positive_f], color="#b7372f", linewidth=1.5, label="Filtered Welch PSD")
    ax.axvline(1.0 / LOWPASS_DAYS, color="0.35", linestyle="--", linewidth=1.0, label="90-day cutoff scale")
    ax.set_title(f"{VARIABLE_LABEL} power spectrum")
    ax.set_xlabel("Frequency (cycles per day)")
    ax.set_ylabel(f"PSD ({UNITS}$^2$ / cycles per day)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(FIGURE_DIR / "assignment1_spectrum_raw_filtered.png", dpi=180)
    plt.close(fig)

    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
