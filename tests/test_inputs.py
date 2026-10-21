"""Reading and cleaning the input: missing values, the 6-REE rule and input checks."""
import numpy as np
import pandas as pd
import pytest

import meltinv
from conftest import REE, write_csv

KEY = ["T_mean", "depth_mean", "basalt_percentage", "misfit_value"]


def test_zero_negative_and_text_values_are_treated_as_missing(in_tmp, grid_sample):
    rows = [grid_sample(k=700)[0], grid_sample(k=1500)[0]]
    blank = pd.DataFrame(rows)
    blank.loc[0, ["La", "Tb"]] = np.nan
    blank.loc[1, "Ho"] = np.nan
    other = blank.astype(object)
    other.loc[0, "La"], other.loc[0, "Tb"], other.loc[1, "Ho"] = 0, -3.0, "bdl"
    blank.to_csv("blank.csv", index=False, float_format="%.17g")
    other.to_csv("other.csv", index=False, float_format="%.17g")
    a = meltinv.invert("blank.csv", group_comparison_figure=False)
    b = meltinv.invert("other.csv", group_comparison_figure=False)
    pd.testing.assert_frame_equal(a[KEY], b[KEY])


def test_samples_with_fewer_than_6_ree_are_skipped(in_tmp, grid_sample, capsys):
    rows = [grid_sample(k=k, location="A")[0] for k in (300, 900, 1600)]
    rows[1] = {**rows[1], **{r: np.nan for r in REE[5:]}}      # second sample keeps only 5 REE
    rows[2] = {**rows[2], **{r: np.nan for r in REE[6:]}}      # third sample keeps exactly 6
    write_csv(rows, "samples.csv")
    out = meltinv.invert("samples.csv", group_comparison_figure=False)
    assert "Skipping A sample 2: 5 usable REE value(s), at least 6 needed." in capsys.readouterr().out
    assert out.sample_number.tolist() == [1, 3]                 # original numbering is kept


def test_no_sample_with_enough_ree_raises(in_tmp, grid_sample):
    row = {**grid_sample()[0], **{r: np.nan for r in REE[4:]}}
    write_csv([row], "few.csv")
    with pytest.raises(ValueError, match="at least 6 usable REE"):
        meltinv.invert("few.csv", group_comparison_figure=False)


@pytest.mark.parametrize("location", [np.nan, "   "])
def test_empty_location_raises(in_tmp, grid_sample, location):
    row = {**grid_sample()[0], "location": location}
    write_csv([row], "bad.csv")
    with pytest.raises(ValueError, match="non-empty location"):
        meltinv.invert("bad.csv", group_comparison_figure=False)


def test_missing_columns_raise(in_tmp, grid_sample):
    row = grid_sample()[0]
    row.pop("Lu")
    write_csv([row], "missing.csv")
    with pytest.raises(ValueError, match="Lu"):
        meltinv.invert("missing.csv", group_comparison_figure=False)


def test_unsupported_file_type_raises(in_tmp):
    open("samples.txt", "w").write("location,La\n")
    with pytest.raises(ValueError, match="Unsupported file format"):
        meltinv.invert("samples.txt", group_comparison_figure=False)


def test_group_inversion_is_not_available_on_main(in_tmp, grid_sample):
    write_csv([grid_sample()[0]], "s.csv")
    with pytest.raises(ValueError, match="group_samples"):
        meltinv.invert("s.csv", per_sample=False)


def test_outputs_are_written_to_the_working_directory(in_tmp, grid_sample):
    write_csv([grid_sample(k=k, location="A")[0] for k in (300, 900, 1600)], "s.csv")
    meltinv.invert("s.csv")                       # group comparison figure is on by default
    assert (in_tmp / "inversion_results/s_sample_inversion_log10.xlsx").exists()
    assert (in_tmp / "inversion_results/s_sample_inversion_summary_log10.xlsx").exists()
    assert (in_tmp / "inversion_figures/group_comparison_A_N3.png").exists()
