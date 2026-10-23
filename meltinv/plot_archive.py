"""Versioned, pickle-free plotting archives; no model grids are needed to read them.

Only the three displayed misfit slices are saved, not the full search cube.
Arrays retain their original precision. The archive is published atomically after
successful inversion; a failed run does not replace an existing archive.
"""
from io import BytesIO
import json
import os
from pathlib import Path
import tempfile
from zipfile import ZipFile, ZIP_DEFLATED

import numpy as np
import pandas as pd

REE = ['La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd',
       'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
SUMMARY_COLUMNS = ['location', 'depth_mean_cali', 'T_mean_cali',
                   'basalt_percentage', 'misfit_value']


def _write_arrays(archive, name, arrays):
    buffer = BytesIO()
    np.savez(buffer, **arrays)
    archive.writestr(name, buffer.getvalue())


def _read_arrays(archive, name):
    with np.load(BytesIO(archive.read(name)), allow_pickle=False) as arrays:
        return {key: arrays[key] for key in arrays.files}


class PlotArchiveWriter:
    """Internal streaming writer. Keep at most one sample's diagnostics in memory."""

    def __init__(self, path, settings):
        self.path = Path(path).expanduser().resolve()
        self.manifest = {'format_version': 1, 'settings': settings, 'samples': []}

    def __enter__(self):
        from . import __version__
        self.manifest['meltinv_version'] = __version__
        self.path.parent.mkdir(parents=True, exist_ok=True)
        fd, self.temp = tempfile.mkstemp(dir=self.path.parent, suffix='.tmp')
        os.close(fd)
        self.archive = ZipFile(self.temp, 'w', compression=ZIP_DEFLATED)
        return self

    def add_sample(self, df, location, count, result):
        cube = result['full_3d_grids']
        p_idx, t_idx = result['min_idx']
        # Store exactly the arrays consumed by plot_results, including uncorrected
        # fit annotations. Group plots use the calibrated per-sample summary below.
        arrays = {key: np.asarray(result[key]) for key in [
            'min_idx', 'misfit_value', 'global_flat_index', 'global_min_value',
            'basalt_idx', 'enrichment_values', 'mean_ree_global', 'mean_ree',
            'low_ree', 'high_ree', 'T_mean', 'depth_mean', 'basalt_percentage',
            'misfit_grid', 'ratio_of_local_global_misfit']}
        arrays.update(grid_shape=np.asarray(cube.shape),
                      temperature_slice=cube[:, p_idx, :],
                      depth_slice=cube[:, :, t_idx],
                      observations=df.loc[df['location'] == location, REE].to_numpy(dtype=float))
        ref = result['litho_thickness_ref']
        if ref is not None:
            arrays['litho_thickness_ref'] = np.asarray(ref)
        name = f'samples/{len(self.manifest["samples"]):06d}.npz'
        _write_arrays(self.archive, name, arrays)
        self.manifest['samples'].append({'file': name, 'location': location, 'count': int(count)})

    def add_summary(self, frame):
        arrays = {key: frame[key].to_numpy(dtype=str if key == 'location' else float)
                  for key in SUMMARY_COLUMNS}
        _write_arrays(self.archive, 'summary.npz', arrays)

    def __exit__(self, exc_type, exc, traceback):
        try:
            if exc_type is None:
                self.archive.writestr('manifest.json', json.dumps(self.manifest, indent=2))
            self.archive.close()
            if exc_type is None:
                os.replace(self.temp, self.path)
        finally:
            self.archive.close()
            if os.path.exists(self.temp):
                os.unlink(self.temp)


def regenerate_figures(archive_path, output_dir='inversion_figures', *,
                       individual_sample_figures=True, group_comparison_figure=True):
    """Render saved diagnostic and group figures without inputs or an inversion.

    Requires a plot_archive produced by invert(). Plot styling remains controlled
    by meltinv.plotting; exact pixel reproduction also needs matching Matplotlib,
    fonts and rendering environment. No pickle objects are loaded.
    """
    from .plotting import plot_results, plot_group_comparison
    with ZipFile(archive_path) as archive:
        manifest = json.loads(archive.read('manifest.json'))
        if manifest['format_version'] != 1:
            raise ValueError('Unsupported meltinv plot archive version')
        if individual_sample_figures:
            for sample in manifest['samples']:
                result = _read_arrays(archive, sample['file'])
                observations = result.pop('observations')
                result = {key: value.item() if value.ndim == 0 else value
                          for key, value in result.items()}
                result.setdefault('litho_thickness_ref', None)
                df = pd.DataFrame(observations, columns=REE)
                df['location'] = sample['location']
                plot_results(df, sample['location'], sample['count'], result, output_dir)
        if group_comparison_figure:
            frame = pd.DataFrame(_read_arrays(archive, 'summary.npz'))
            plot_group_comparison(frame, output_dir=output_dir)
