"""Run from the repository root: PYTHONPATH=. python review/reproduce_audit.py."""
import ast
import json
import subprocess
import numpy as np
import pandas as pd
from meltinv.inversion import load_all_grids, source_composition_mix

BASE = '48fc8e473acb13cfec6722468a3cecc5750b2b1c'
source = subprocess.check_output(['git', 'show', f'{BASE}:meltinv/inversion.py'], text=True)
# Execute only the baseline numerical functions, avoiding its eager meltPT import.
names = {'compute_total_ree_misfit', 'source_composition_mix', 'apply_mixing'}
nodes = [n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name in names]
namespace = {'np': np, 'pd': pd}
exec(compile(ast.Module(body=nodes, type_ignores=[]), '<baseline>', 'exec'), namespace)
g = load_all_grids()
valid = np.argwhere(np.isfinite(g['simulation']).all(axis=0))
r, c = valid[500]
observed = g['simulation'][:, r, c]
df = pd.DataFrame([observed], columns=g['ree_variables'])
old = namespace['compute_total_ree_misfit'](df, [True], 0, g, 1)
result = {'baseline_commit': BASE, 'grid_shape': list(g['simulation'].shape[1:]),
          'valid_cells': len(valid), 'exact_input_cell': [int(r), int(c)],
          'baseline_selected_cell': [int(i) for i in old[2]],
          'baseline_selected_misfit': float(old[0][old[2]]),
          'baseline_misfit_at_exact_cell': float(old[0][r,c])}
# Illustrative same-grid perturbations, not validation against natural samples.
rng = np.random.default_rng(20261019)
scales = np.arange(0, 15.1, .1)
models = np.stack([g['simulation'] * np.array(list(source_composition_mix(b).values()))[:,None,None]
                   for b in scales])
logs = np.log10(models)
comparisons = []
for idx in rng.choice(len(valid), 24, replace=False):
    r, c = valid[idx]
    observation = np.log10(g['simulation'][:, r, c]) + rng.normal(0, .03, 14)
    residual = logs - observation[None,:,None,None]
    l1 = np.mean(abs(residual), axis=1)
    l2 = np.sqrt(np.mean(residual**2, axis=1))
    a = np.unravel_index(np.nanargmin(l1), l1.shape)
    b = np.unravel_index(np.nanargmin(l2), l2.shape)
    comparisons.append({'mae_index': list(map(int,a)), 'rms_index': list(map(int,b)),
                        'temperature_difference_K': int((b[2]-a[2])*5),
                        'thickness_difference_km': int((b[1]-a[1])*2),
                        'enrichment_difference_percentage_points': float(scales[b[0]]-scales[a[0]])})
envelope_findings = []
for i, element in enumerate(g['ree_variables']):
    mean, low, high = g['simulation'][i], g['lower'][i], g['upper'][i]
    violation = np.isfinite(mean) & ((low > mean) | (mean > high))
    if violation.any():
        relative = np.maximum(low-mean, mean-high)/mean
        row, col = np.unravel_index(np.nanargmax(relative), relative.shape)
        envelope_findings.append({'element': element, 'count': int(violation.sum()),
            'max_relative_violation': float(relative[row,col]),
            'worst_cell': [int(row),int(col)], 'min_mean_max': list(map(float,
            [low[row,col],mean[row,col],high[row,col]]))})
result['envelope_violations'] = envelope_findings
result['synthetic_noise_dex'] = .03
result['synthetic_seed'] = 20261019
result['synthetic_changed_optima'] = sum(x['mae_index'] != x['rms_index'] for x in comparisons)
result['synthetic_comparisons'] = comparisons
with open('review/audit_results.json','w') as out:
    json.dump(result,out,indent=2)
print(json.dumps({k:v for k,v in result.items() if k != 'synthetic_comparisons'}, indent=2))
