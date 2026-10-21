"""
Reference-result (regression) test on the bundled example.

The whole inversion is run on the olivine-corrected examples (tests/reference/examples_corrected.csv)
and every result is compared with the stored reference. If this test fails after a change you
INTENDED to make, regenerate the reference with `python tests/update_reference.py` and explain why
in the commit message. If you did not intend to change the results, the failure points to a bug.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import meltinv
from update_reference import SAMPLE_COLS, SUMMARY_COLS, corrected_input

REF = Path(__file__).resolve().parent / "reference"
GRID_COLS = ["T_mean", "depth_mean", "basalt_percentage"]     # grid values: must match exactly


@pytest.fixture(scope="module")
def results(tmp_path_factory, monkeypatch_module):
    workdir = tmp_path_factory.mktemp("reference")
    monkeypatch_module.chdir(workdir)
    samples = meltinv.invert(str(REF / "examples_corrected.csv"), group_comparison_figure=False)
    summary = pd.read_excel(workdir / "inversion_results/examples_sample_inversion_summary_log10.xlsx")
    return samples, summary


@pytest.fixture(scope="module")
def monkeypatch_module():
    mp = pytest.MonkeyPatch()
    yield mp
    mp.undo()


def test_per_sample_results_match_reference(results):
    samples, _ = results
    ref = pd.read_csv(REF / "examples_samples.csv")
    got = samples[SAMPLE_COLS].reset_index(drop=True)
    assert len(got) == len(ref)
    assert got.location.tolist() == ref.location.tolist()
    assert got.sample_number.tolist() == ref.sample_number.tolist()
    for c in GRID_COLS:
        np.testing.assert_allclose(got[c], ref[c], rtol=0, atol=1e-9, err_msg=c)
    for c in set(SAMPLE_COLS) - set(GRID_COLS) - {"location", "sample_number"}:
        np.testing.assert_allclose(got[c], ref[c], rtol=1e-6, err_msg=c)


def test_location_summary_matches_reference(results):
    _, summary = results
    ref = pd.read_csv(REF / "examples_summary.csv")
    got = summary[SUMMARY_COLS].reset_index(drop=True)
    assert got.location.tolist() == ref.location.tolist()
    assert got["count"].tolist() == ref["count"].tolist()
    for c in set(SUMMARY_COLS) - {"location", "count"}:
        np.testing.assert_allclose(got[c], ref[c], rtol=1e-6, err_msg=c)


@pytest.mark.slow
def test_olivine_correction_matches_reference(tmp_path, monkeypatch):
    """Runs meltPT on examples.xlsx (slow, and needs meltPT). Run with: pytest -m slow"""
    pytest.importorskip("meltPT")
    monkeypatch.chdir(tmp_path)                    # restored after the test
    got = corrected_input(tmp_path).reset_index(drop=True)
    ref = pd.read_csv(REF / "examples_corrected.csv")
    assert got.correction_success.tolist() == ref.correction_success.tolist()
    np.testing.assert_allclose(got.ol_added, ref.ol_added, rtol=1e-6)
    for c in ["La", "Sm", "Yb"]:
        np.testing.assert_allclose(got[c], ref[c], rtol=1e-6, err_msg=c)
