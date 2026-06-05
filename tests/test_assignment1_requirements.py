import numpy as np
import pytest

from spectra_filtering.analysis import summary_stats
from spectra_filtering.spectra import parseval_ratio, welch_psd


def test_summary_stats_mean_ignores_missing_values():
    stats = summary_stats(np.array([1.0, np.nan, 3.0, 5.0]))

    assert stats["mean"] == 3.0
    assert stats["n_missing"] == 1


def test_welch_psd_satisfies_parseval_for_regular_series():
    rng = np.random.default_rng(42)
    x = rng.normal(size=4096)

    freq, psd = welch_psd(x, dt_days=0.5, segment_length=1024)

    assert parseval_ratio(x, freq, psd) == pytest.approx(1.0, rel=0.08)
