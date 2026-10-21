"""
Regenerate the reference results used by tests/test_reference.py.

Only run this after an INTENDED change to the results (e.g. a new misfit or a new
calibration), and say why in the commit message. Run from the repository root:

    python tests/update_reference.py           # re-invert the stored corrected input
    python tests/update_reference.py --correct # also redo the olivine correction (needs meltPT)
"""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
import pandas as pd

import meltinv

REF = Path(__file__).resolve().parent / "reference"
REE = ['La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
SAMPLE_COLS = ["location", "sample_number", "T_mean", "depth_mean", "basalt_percentage", "T_mean_cali",
               "depth_mean_cali", "pres_mean_cali", "F_mean", "M_mean", "misfit_value"]
SUMMARY_COLS = ["location", "count", "T_mean_cali", "T_mean_cali_err_minus", "T_mean_cali_err_plus",
                "depth_mean_cali", "basalt_percentage", "misfit_value"]
CORRECTED_COLS = ["location", "Sample", "thickness", *REE, "ol_added", "correction_success"]


def corrected_input(workdir):
    """Olivine-correct examples.xlsx with meltPT and keep only the columns the inversion needs."""
    examples = Path(__file__).resolve().parents[1] / "examples.xlsx"
    os.chdir(workdir)
    pd.read_excel(examples).to_excel("examples.xlsx", index=False)
    meltinv.correction("examples.xlsx")
    return pd.read_excel("inversion_results/examples_corrected.xlsx")[CORRECTED_COLS]


def run_reference(corrected_csv, workdir):
    os.chdir(workdir)
    samples = meltinv.invert(str(corrected_csv), group_comparison_figure=False)
    # meltinv drops "_corrected" from the input name when naming its outputs
    summary = pd.read_excel("inversion_results/examples_sample_inversion_summary_log10.xlsx")
    return samples[SAMPLE_COLS], summary[SUMMARY_COLS]


if __name__ == "__main__":
    REF.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        if "--correct" in sys.argv:
            corrected_input(tmp).to_csv(REF / "examples_corrected.csv", index=False, float_format="%.17g")
            print("updated", REF / "examples_corrected.csv")
        samples, summary = run_reference(REF / "examples_corrected.csv", tmp)
    samples.to_csv(REF / "examples_samples.csv", index=False, float_format="%.17g")
    summary.to_csv(REF / "examples_summary.csv", index=False, float_format="%.17g")
    print("updated", REF / "examples_samples.csv")
    print("updated", REF / "examples_summary.csv")
