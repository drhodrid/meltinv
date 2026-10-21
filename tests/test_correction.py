"""Olivine fractionation correction (the part that does not need meltPT itself)."""
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from meltinv.frac_correction import save_fraction
from conftest import REE, MAJORS


def primary_output(ol_added, complete):
    """Mimic meltPT's Suite.primary: olivine added and primary major elements per sample."""
    return pd.DataFrame({"ol_added": ol_added,
                         **{f"{m}_primary_wt": [1.0 if ok else np.nan for ok in complete] for m in MAJORS}})


def test_successful_correction_dilutes_ree_by_added_olivine(in_tmp):
    sample = pd.DataFrame({c: [10.0, 10.0] for c in ["location", *REE, *MAJORS]})
    sample["location"] = "A"
    out = save_fraction("s.csv", sample, primary_output([0.2, 0.0], [True, True]))
    assert out.La.tolist() == pytest.approx([8.0, 10.0])
    assert out.correction_success.tolist() == [True, True]
    assert "SiO2_primary" in out.columns
    assert (in_tmp / "inversion_results/s_corrected.xlsx").exists()


def test_failed_correction_gives_missing_ree_not_zero_olivine(in_tmp):
    # Regression (PR #1): failed corrections used to be treated as 0 % olivine added,
    # so uncorrected REE entered the inversion without any warning.
    sample = pd.DataFrame({c: [10.0, 10.0, 10.0] for c in ["location", *REE, *MAJORS]})
    sample["location"] = "A"
    primary = primary_output([0.2, np.nan, 0.3], [True, False, False])
    with pytest.warns(RuntimeWarning, match="failed for 2 sample"):
        out = save_fraction("s.csv", sample, primary)
    assert out.La.iloc[0] == pytest.approx(8.0)
    assert out.loc[1:, REE].isna().all().all()
    assert out.correction_success.tolist() == [True, False, False]
    assert sample.La.eq(10.0).all()                 # the input table is not modified


def test_olivine_addition_outside_0_to_1_counts_as_failed(in_tmp):
    sample = pd.DataFrame({c: [10.0, 10.0] for c in ["location", *REE, *MAJORS]})
    sample["location"] = "A"
    with pytest.warns(RuntimeWarning):
        out = save_fraction("s.csv", sample, primary_output([1.0, -0.1], [True, True]))
    assert out.correction_success.tolist() == [False, False]


def test_importing_meltinv_does_not_load_meltpt():
    # meltPT is only needed when a correction is requested.
    code = "import sys, meltinv; print('meltPT' in sys.modules)"
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert result.stdout.strip() == "False"
