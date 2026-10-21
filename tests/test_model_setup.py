"""Packaged model grids, source compositions, enrichment ranges and the T-P correction."""
import numpy as np
import pytest

from meltinv.inversion import calibration, get_enrichment_values, source_composition_mix
from conftest import REE


def test_packaged_grids_have_the_expected_shape(grids):
    assert grids["ree_variables"] == REE
    assert grids["simulation"].shape == (14, 53, 58)        # REE x thickness x temperature
    for key in ["lower", "upper"]:
        assert grids[key].shape == grids["simulation"].shape
    for key in ["pres_mean", "F_mean", "M_mean"]:
        assert grids[key].shape == (53, 58)


def test_enrichment_search_ranges():
    default = get_enrichment_values("Hawaii", None)
    assert len(default) == 151 and default[0] == 0 and default[-1] == pytest.approx(15.0)
    assert len(get_enrichment_values("Iceland (NNVZ)", ["Iceland"])) == 171     # adds -100 to -5
    assert len(get_enrichment_values("Hawaii", ["Iceland"])) == 151
    assert len(get_enrichment_values("Hawaii", "all")) == 171


def test_source_mixing_end_members():
    # Primitive mantle (McDonough & Sun 1995), basalt (Sun & McDonough 1989), depleted mantle (Salters & Stracke 2004)
    assert all(v == pytest.approx(1.0) for v in source_composition_mix(0).values())
    assert source_composition_mix(-100)["La"] == pytest.approx(0.234 / 0.648)
    assert source_composition_mix(15)["Lu"] == pytest.approx((0.15 * 0.455 + 0.85 * 0.0675) / 0.0675)


def test_t_p_correction_is_zero_without_enrichment_and_grows_with_it():
    T, P, h = np.array([150.0] * 3), np.array([100.0] * 3), np.array([60.0] * 3)
    T_cali, h_cali = calibration(T, P, h, np.array([0.0, 5.0, 10.0]))
    assert T_cali[0] == pytest.approx(150.0) and h_cali[0] == pytest.approx(60.0)
    assert T_cali[2] < T_cali[1] < 150.0          # cooler with more pyroxenite
    assert h_cali[2] > h_cali[1] > 60.0           # thicker with more pyroxenite


def test_depleted_sources_get_no_t_p_correction():
    T_cali, h_cali = calibration(np.array([150.0]), np.array([100.0]), np.array([60.0]), np.array([-50.0]))
    assert T_cali[0] == pytest.approx(150.0) and h_cali[0] == pytest.approx(60.0)
