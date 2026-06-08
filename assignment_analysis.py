"""Run Assignment 1 analysis for a NAC AMOC time series."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from amocatlas import read
from scipy import signal
from scipy.stats import chi2

from spectra_filtering.analysis import summary_stats
from spectra_filtering.data_io import fill_gaps
from spectra_filtering.filters import tukey_lowpass
from spectra_filtering.spectra import parseval_ratio, welch_psd


FIGURE_DIR = Path("figures")
OUTPUT_DIR = Path("outputs")

ARRAY_NAME = "NAC"
VARIABLE = "TRANS_NAC_PROXY"
VARIABLE_LABEL = "North Atlantic Current transport proxy"
UNITS = "Sv"
SEGMENT_LENGTH = 64
OVERLAP = 0.0
LOWPASS_DAYS = 5 * 365.25
CONFIDENCE_LEVEL = 0.95
DAYS_PER_YEAR = 365.25


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


def _welch_segments(n: int, segment_length: int, overlap: float) -> int:
    step = int(round(segment_length * (1.0 - overlap)))
    if step <= 0:
        raise ValueError("overlap leaves no positive Welch step")
    return 1 + max(0, (n - segment_length) // step)


def _confidence_interval_factors(dof: float, level: float) -> tuple[float, float]:
    alpha = 1.0 - level
    lower = dof / chi2.ppf(1.0 - alpha / 2.0, dof)
    upper = dof / chi2.ppf(alpha / 2.0, dof)
    return float(lower), float(upper)


def _window_response(window: str, length: int, dt_days: float) -> tuple[np.ndarray, np.ndarray]:
    if window == "tukey":
        weights = signal.windows.tukey(length, alpha=0.5)
    elif window == "boxcar":
        weights = np.ones(length)
    else:
        raise ValueError(f"unsupported window: {window}")
    weights = weights / np.sum(weights)
    n_fft = 8192
    response = np.abs(np.fft.rfft(weights, n=n_fft)) ** 2
    response = response / response[0]
    freq = np.fft.rfftfreq(n_fft, d=dt_days)
    return freq, response


def _load_series() -> tuple[np.ndarray, np.ndarray, float]:
    ds = read.nac()
    da = ds[VARIABLE]
    time = da["TIME"].values
    values = np.squeeze(da.values.astype("float64"))
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
    n_segments = _welch_segments(len(filled), SEGMENT_LENGTH, OVERLAP)
    dof = 2 * n_segments
    ci_lower_factor, ci_upper_factor = _confidence_interval_factors(dof, CONFIDENCE_LEVEL)

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
        "welch_segments": n_segments,
        "welch_degrees_of_freedom": dof,
        "confidence_level": CONFIDENCE_LEVEL,
        "confidence_lower_factor": ci_lower_factor,
        "confidence_upper_factor": ci_upper_factor,
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
    ax.plot(time, filtered, color="#b7372f", linewidth=1.8, label="5-year Tukey low-pass")
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
    freq_year = freq * DAYS_PER_YEAR
    freq_year_f = freq_f * DAYS_PER_YEAR
    psd_year = psd / DAYS_PER_YEAR
    psd_year_f = psd_f / DAYS_PER_YEAR
    ax.loglog(freq_year[positive], psd_year[positive], color="#31688e", linewidth=1.1, label="Raw Welch PSD")
    ax.fill_between(
        freq_year[positive],
        psd_year[positive] * ci_lower_factor,
        psd_year[positive] * ci_upper_factor,
        color="#31688e",
        alpha=0.18,
        linewidth=0,
        label="95% chi-squared CI",
    )
    ax.loglog(freq_year_f[positive_f], psd_year_f[positive_f], color="#b7372f", linewidth=1.5, label="Filtered Welch PSD")
    ax.axvline(DAYS_PER_YEAR / LOWPASS_DAYS, color="0.35", linestyle="--", linewidth=1.0, label="5-year cutoff scale")
    ax.set_title(f"{VARIABLE_LABEL} power spectrum")
    ax.set_xlabel("Frequency (cycles per year)")
    ax.set_ylabel(f"PSD ({UNITS}$^2$ / cycles per year)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(FIGURE_DIR / "assignment1_spectrum_raw_filtered.png", dpi=180)
    plt.close(fig)

    freq_resp_t, response_t = _window_response("tukey", lp_window, dt_days)
    freq_resp_b, response_b = _window_response("boxcar", lp_window, dt_days)
    freq_resp_year_t = freq_resp_t * DAYS_PER_YEAR
    freq_resp_year_b = freq_resp_b * DAYS_PER_YEAR
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    valid_t = freq_resp_t > 0
    valid_b = freq_resp_b > 0
    ax.semilogx(freq_resp_year_t[valid_t], response_t[valid_t], color="#b7372f", linewidth=1.6, label="Tukey window")
    ax.semilogx(freq_resp_year_b[valid_b], response_b[valid_b], color="#31688e", linewidth=1.2, linestyle="--", label="Boxcar window")
    ax.axvline(DAYS_PER_YEAR / LOWPASS_DAYS, color="0.35", linestyle="--", linewidth=1.0, label="5-year scale")
    ax.set_ylim(-0.03, 1.05)
    ax.set_title("Low-pass window power response")
    ax.set_xlabel("Frequency (cycles per year)")
    ax.set_ylabel("Normalised power response")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(FIGURE_DIR / "assignment1_filter_response.png", dpi=180)
    plt.close(fig)

    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
