"""A saved run must reproduce both plot types without touching the inversion."""
import subprocess
import sys
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pytest
from meltinv import invert, regenerate_figures
from meltinv.plot_archive import PlotArchiveWriter
from conftest import write_csv


def test_inversion_import_does_not_import_matplotlib():
    subprocess.run([sys.executable, '-c',
                    "import meltinv; import sys; assert 'matplotlib.pyplot' not in sys.modules"],
                   check=True)


def test_archive_reproduces_figures_and_preserves_results(in_tmp, grid_sample, monkeypatch):
    rows = [grid_sample(k=500+i, enrichment=-35, location='A')[0] for i in range(3)]
    rows[0]['Lu'] = np.nan  # missing REE must survive saving and plotting
    write_csv(rows, 'input.csv')
    plain = invert('input.csv', depleted_locations='all', group_comparison_figure=False)
    saved = invert('input.csv', depleted_locations='all', individual_sample_figures=True,
                   plot_archive='run.zip')
    pd.testing.assert_frame_equal(plain, saved, check_exact=True)
    original = {p.name: p.read_bytes() for p in (in_tmp/'inversion_figures').glob('*.png')}
    assert len(original) == 4
    with ZipFile('run.zip') as z:
        assert len(z.namelist()) == 5
    (in_tmp/'input.csv').unlink()
    def forbidden(*args, **kwargs):
        raise AssertionError('Replot must not load model grids or rerun inversion')
    monkeypatch.setattr('meltinv.inversion.load_all_grids', forbidden)
    monkeypatch.setattr('meltinv.inversion.invert_single_group', forbidden)
    regenerate_figures('run.zip', 'replotted')
    regenerated = {p.name: p.read_bytes() for p in (in_tmp/'replotted').glob('*.png')}
    assert original == regenerated  # same rendering environment, exact PNG bytes


def test_failed_archive_does_not_replace_previous_run(tmp_path):
    path = tmp_path/'run.zip'
    path.write_bytes(b'previous run')
    with pytest.raises(RuntimeError):
        with PlotArchiveWriter(path, {}):
            raise RuntimeError('failed inversion')
    assert path.read_bytes() == b'previous run'
    assert list(tmp_path.iterdir()) == [path]


def test_archive_can_be_saved_without_rendering(in_tmp, grid_sample, monkeypatch):
    write_csv([grid_sample()[0]], 'input.csv')
    def forbidden(*args, **kwargs):
        raise AssertionError('Saving an archive must not render figures')
    monkeypatch.setattr('meltinv.plotting.plot_results', forbidden)
    monkeypatch.setattr('meltinv.plotting.plot_group_comparison', forbidden)
    invert('input.csv', group_comparison_figure=False, plot_archive='run.zip')
    assert (in_tmp/'run.zip').exists()
    assert not (in_tmp/'inversion_figures').exists()


def test_cli_plot_dispatch(monkeypatch):
    from meltinv.cli import main
    calls = []
    monkeypatch.setattr('meltinv.plot_archive.regenerate_figures',
                        lambda *args, **kwargs: calls.append((args, kwargs)))
    monkeypatch.setattr(sys, 'argv', ['meltinv', 'plot', 'run.zip',
                                     '--output_dir', 'figures', '--no_individual'])
    main()
    assert calls == [(('run.zip', 'figures'), {
        'individual_sample_figures': False, 'group_comparison_figure': True})]
