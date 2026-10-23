"""Location summary: all samples, weighted by 1/misfit."""
import numpy as np
import pandas as pd
import pytest

from meltinv.inversion import aggregate_per_sample_results

COLUMNS = ['T_mean', 'T_mean_cali', 'depth_mean', 'depth_mean_cali', 'pres_mean_cali', 'pres_min_cali',
           'pres_max_cali', 'F_mean', 'F_min', 'F_max', 'M_mean', 'M_min', 'M_max', 'basalt_percentage',
           'misfit_value', 'global_min_value', 'ratio_of_local_global_misfit', 'litho_thickness_ref']


def per_sample(location, T, misfit, enrichment=None):
    n = len(T)
    df = pd.DataFrame({c: np.ones(n) for c in COLUMNS})
    df["location"] = location
    df["T_mean_cali"] = T
    df["misfit_value"] = misfit
    if enrichment is not None:
        df["basalt_percentage"] = enrichment
    return df


def test_all_samples_are_averaged_with_1_over_misfit_weights():
    # Ten samples, so the old "best half" rule would have applied.
    T = np.arange(100.0, 200.0, 10.0)
    misfit = np.linspace(0.01, 0.10, 10)
    out = aggregate_per_sample_results(per_sample("A", T, misfit)).iloc[0]
    w = 1 / (misfit + 1e-6)
    mean = np.average(T, weights=w)
    std = np.sqrt(np.average((T - mean) ** 2, weights=w))
    assert out.T_mean_cali == pytest.approx(mean)
    assert out.T_mean_cali_err_minus == pytest.approx(mean - std)
    assert out.T_mean_cali_err_plus == pytest.approx(mean + std)
    assert out["count"] == 10
    assert "count_used" not in out.index           # removed: every sample is used


def test_each_location_is_summarised_separately():
    df = pd.concat([per_sample("A", [100.0, 120.0], [0.02, 0.02]), per_sample("B", [200.0], [0.03])])
    out = aggregate_per_sample_results(df).set_index("location")
    assert out.loc["A", "T_mean_cali"] == pytest.approx(110.0)
    assert out.loc["B", "T_mean_cali"] == pytest.approx(200.0)
    assert out.loc["B", "T_mean_cali_err_minus"] == pytest.approx(200.0)    # one sample: no spread


def test_depleted_enrichment_bounds_are_in_the_same_units_as_the_mean():
    # Negative (depleted) values stay in percent; only figures compress them.
    out = aggregate_per_sample_results(per_sample("A", [100.0, 100.0], [0.02, 0.02], [-100.0, -50.0])).iloc[0]
    assert out.basalt_percentage == pytest.approx(-75.0)
    assert out.basalt_percentage_err_minus == pytest.approx(-100.0)
    assert out.basalt_percentage_err_plus == pytest.approx(-50.0)
