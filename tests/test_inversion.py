import numpy as np
import pandas as pd
import pytest
from meltinv import invert
from meltinv.inversion import (compute_total_ree_misfit, load_all_grids,
    source_composition_mix, get_enrichment_values, calibration,
    aggregate_per_sample_results)
from meltinv.frac_correction import save_fraction


def objective(models, observations, metric='mae'):
    names = ['La', 'Ce']
    df = pd.DataFrame([observations], columns=names)
    grids = {'ree_variables': names, 'simulation': np.asarray(models, dtype=float)}
    return compute_total_ree_misfit(df, [True], 0, grids, 1, metric)


def test_objectives_can_select_different_models():
    # Residuals A=(0, 2), B=(1.1, 1.1): L1 chooses A, L2 chooses B.
    models = np.power(10., [[[0, 1.1]], [[2, 1.1]]])
    mae = objective(models, [1, 1])
    paper = objective(models, [1, 1], 'paper_l2')
    rms = objective(models, [1, 1], 'rms')
    assert mae[2] == (0, 0)
    assert paper[2] == rms[2] == (0, 1)
    np.testing.assert_allclose(paper[0], [[1, 1.1 / np.sqrt(2)]])
    np.testing.assert_allclose(rms[0], [[np.sqrt(2), 1.1]])


def test_exact_fit_and_constant_grid_are_valid():
    assert objective([[[1, 2]], [[1, 2]]], [1, 1])[2] == (0, 0)
    assert objective([[[2, 2]], [[2, 2]]], [1, 1])[2] == (0, 0)


@pytest.mark.parametrize('bad', [0, -1, np.inf, np.nan])
def test_invalid_second_element_excludes_model_cell(bad):
    result = objective([[[1, 2]], [[bad, 2]]], [1, 1])
    assert result[2] == (0, 1)
    assert np.isnan(result[0][0, 0])


def test_source_ranges_and_table_values():
    assert len(get_enrichment_values('A', None)) == 151
    assert len(get_enrichment_values('A', 'all')) == 171
    assert source_composition_mix(0)['La'] == 1
    assert source_composition_mix(-100)['La'] == pytest.approx(.234/.648)
    assert source_composition_mix(15)['Lu'] == pytest.approx((.15*.455+.85*.0675)/.0675)


def test_correction_equations_and_depleted_identity():
    t, h = calibration([200, 200], [99, 99], np.array([60, 60]), [-50, 10])
    x = .9 / 1.8
    offset = 920 + 120*3 - (1101.3 + 141.4*3 - 6.8*9)
    np.testing.assert_allclose(t, [200, 200 + .5*x*offset])
    np.testing.assert_allclose(h, [60, 60 - 33*.5*x*offset/120])


def test_summary_enrichment_bounds_use_same_units_as_mean():
    columns = ['T_mean', 'T_mean_cali', 'depth_mean', 'depth_mean_cali',
        'pres_mean_cali', 'pres_min_cali', 'pres_max_cali', 'F_mean', 'F_min',
        'F_max', 'M_mean', 'M_min', 'M_max', 'basalt_percentage', 'misfit_value',
        'global_min_value', 'ratio_of_local_global_misfit', 'litho_thickness_ref']
    df = pd.DataFrame({c: [1., 1.] for c in columns})
    df['location'] = 'A'; df['basalt_percentage'] = [-100., -50.]
    out = aggregate_per_sample_results(df).iloc[0]
    assert out.basalt_percentage == -75
    assert out.basalt_percentage_err_minus == -100
    assert out.basalt_percentage_err_plus == -50
    assert out['count_used'] == 2


def test_failed_fractionation_is_not_zero_addition(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    ree = load_all_grids()['ree_variables']
    majors = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'CaO', 'MgO', 'MnO', 'K2O', 'Na2O', 'P2O5']
    sample = pd.DataFrame({c: [10., 10., 10.] for c in ree + majors})
    primary = pd.DataFrame({c+'_primary_wt': [1., np.nan, np.nan] for c in majors})
    primary['ol_added'] = [.2, np.nan, .51]
    with pytest.warns(RuntimeWarning, match='failed for 2'):
        out = save_fraction('test.csv', sample, primary)
    assert out.La.iloc[0] == 8
    assert out.La.iloc[1:].isna().all()
    assert out.correction_success.tolist() == [True, False, False]
    assert sample.La.eq(10).all()


@pytest.mark.parametrize('metric', ['mae', 'paper_l2', 'rms'])
def test_packaged_grid_exact_recovery_end_to_end(tmp_path, monkeypatch, metric):
    monkeypatch.chdir(tmp_path)
    grids = load_all_grids()
    row, col = np.argwhere(np.isfinite(grids['simulation']).all(axis=0))[500]
    data = {r: [grids['simulation'][i, row, col]] for i, r in enumerate(grids['ree_variables'])}
    data['location'] = ['synthetic']
    pd.DataFrame(data).to_csv('synthetic.csv', index=False, float_format='%.17g')
    out = invert('synthetic.csv', group_comparison_figure=False, misfit_metric=metric)
    assert out.T_mean.iloc[0] == 65 + 5*col
    assert out.depth_mean.iloc[0] == 16 + 2*row  # existing convention, pending grid provenance
    assert out.basalt_percentage.iloc[0] == 0
    assert out.misfit_value.iloc[0] < 1e-14
    assert out.misfit_metric.iloc[0] == metric
    assert (tmp_path/'inversion_results/synthetic_sample_inversion_summary_log10.xlsx').exists()


def test_infinite_observations_are_missing():
    result = objective([[[1, 2]], [[1, 2]]], [np.inf, 1])
    assert result[2] == (0, 0)
    assert result[0][0, 0] == 0


def test_no_model_support_returns_no_minimum():
    result = objective([[[np.nan]], [[1]]], [1, 1])
    assert result[2] is None


def test_bad_metric_is_rejected():
    with pytest.raises(ValueError, match='misfit_metric'):
        objective([[[1]], [[1]]], [1, 1], 'unknown')


def test_legacy_group_preserves_upstream_per_element_rms():
    df = pd.DataFrame({'La': [1., 100.]})
    grids = {'ree_variables': ['La'], 'simulation': np.array([[[1., 10.]]])}
    result = compute_total_ree_misfit(df, [True, True], 0, grids, 2)
    np.testing.assert_allclose(result[0], [[np.sqrt(2), 1.]])
