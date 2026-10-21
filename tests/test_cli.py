"""Command-line interface."""
import subprocess
import sys

from conftest import write_csv


def run_cli(*args, cwd=None):
    return subprocess.run([sys.executable, "-m", "meltinv.cli", *args], cwd=cwd,
                          capture_output=True, text=True)


def test_help_lists_both_commands():
    result = run_cli("--help")
    assert result.returncode == 0
    assert "correction" in result.stdout and "invert" in result.stdout


def test_invert_runs_from_the_command_line(in_tmp, grid_sample):
    write_csv([grid_sample(k=k, location="A")[0] for k in (300, 900, 1600)], "s.csv")
    result = run_cli("invert", "s.csv", "--no_group_comparison_figure", cwd=in_tmp)
    assert result.returncode == 0, result.stderr
    assert (in_tmp / "inversion_results/s_sample_inversion_log10.xlsx").exists()
    assert not (in_tmp / "inversion_figures").exists()


def test_old_flags_are_still_accepted(in_tmp, grid_sample):
    write_csv([grid_sample(k=300)[0]], "s.csv")
    result = run_cli("invert", "s.csv", "--per_sample", "--make_figures", "--no_group_comparison_figure",
                     cwd=in_tmp)
    assert result.returncode == 0, result.stderr
    assert (in_tmp / "inversion_figures/ree_comparison_synthetic__sample1_log10_N1.png").exists()
