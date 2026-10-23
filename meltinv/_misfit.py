"""Prepared log-concentrations for the exhaustive search (no scientific settings).

Model arrays have axes (REE, thickness row, temperature column). The cache is
local to one inversion run, so it cannot reuse predictions from another grid.
"""
import numpy as np
import pandas as pd


class ModelLogCache:
    """Cache each enrichment's log predictions and positive/finite mask.

    Keep log10(concentration * factor), rather than adding separate logarithms:
    this preserves the original floating-point operations, including exact fits.
    At 171 enrichments, 14 REEs and 53 x 58 cells the arrays occupy about 63 MiB.
    Treat the source grids as immutable for the lifetime of this cache.
    """

    def __init__(self, grids, mixing_function):
        self.grids = grids
        self.mixing_function = mixing_function
        self._entries = {}

    def get(self, enrichment):
        if enrichment not in self._entries:
            factors = self.mixing_function(enrichment)
            scaled = np.asarray([
                np.asarray(grid, dtype=float) * float(factors[ree])
                for grid, ree in zip(self.grids['simulation'], self.grids['ree_variables'])
            ])
            valid = np.isfinite(scaled) & (scaled > 0)
            logs = np.full_like(scaled, np.nan)
            logs[valid] = np.log10(scaled[valid])
            self._entries[enrichment] = logs, valid, factors
        return self._entries[enrichment]


def prepare_observations(df, mask, ree_variables, min_sample_count):
    """Extract usable log10 ppm values once, retaining element and sample order."""
    observations = []
    for i, ree in enumerate(ree_variables):
        values = pd.to_numeric(df.loc[mask, ree], errors='coerce').to_numpy(dtype=float)
        values = values[np.isfinite(values) & (values > 0)]
        if len(values) >= min_sample_count:
            observations.append((i, np.log10(values)))
    return observations


def evaluate_misfit(observations, model_logs, model_valid, *, diagnostic_threshold=True):
    """Return the N-3 log10 residual grid, optional percentile and first minimum.

    Any invalid prediction for a measured element excludes that model cell.
    Summation order matches the original element-then-sample loop. Ties select
    the first cell in row-major order; zero residuals remain valid minima.
    """
    n = sum(len(values) for _, values in observations)
    if not observations or n <= 3:
        return None, None, None
    valid = np.ones(model_logs.shape[1:], dtype=bool)
    squared = np.zeros(model_logs.shape[1:], dtype=float)
    for i, values in observations:
        usable = model_valid[i]
        valid &= usable
        for value in values:
            squared[usable] += (model_logs[i][usable] - value) ** 2
    misfit = np.where(valid, np.sqrt(squared / (n - 3)), np.nan)
    if not valid.any():
        return misfit, None, None
    threshold = np.percentile(misfit[valid], 10) if diagnostic_threshold else None
    index = np.unravel_index(np.nanargmin(misfit), misfit.shape)
    return misfit, threshold, index
