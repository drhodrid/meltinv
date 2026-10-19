# meltinv

**meltinv** is a Python package that uses the Rare Earth Element (REE) compositions of primary magmas to jointly invert for:

- **Mantle excess temperature** (mantle potential temperature above ambient, in K)
- **Lithospheric thickness** (km), and the depth range over which melting occurs
- **Relative source enrichment**, meaning how much the mantle source is enriched (or depleted) relative to primitive mantle

It also reports the melt fraction (`F`) and melting rate (`M`) of the best-fitting model.

The package can be used from the **command line** or from **Python**. This README is written as a step-by-step tutorial. Work through it once with the bundled example file, then swap in your own data.

---

## Contents

1. [Installation](#1-installation)
2. [Tutorial: your first inversion](#2-tutorial-your-first-inversion)
   - [Step 1: Prepare the input file](#step-1-prepare-the-input-file)
   - [Step 2: (Optional) Fractionation correction](#step-2-optional-fractionation-correction)
   - [Step 3: Run the inversion](#step-3-run-the-inversion)
   - [Step 4: Reading the results](#step-4-reading-the-results)
   - [Step 5: Reading the figures](#step-5-reading-the-figures)
3. [Going further](#3-going-further)
   - [Depleted mantle sources](#depleted-mantle-sources)
   - [Inverting each sample individually](#inverting-each-sample-individually)
   - [Using the results in Python](#using-the-results-in-python)
4. [Reference](#4-reference)
5. [How the inversion works](#5-how-the-inversion-works)
6. [Troubleshooting](#6-troubleshooting)
7. [Citation](#7-citation)

---

## 1. Installation

meltinv requires **Python ≥ 3.11**. We recommend installing it in a fresh virtual environment:

```bash
git clone https://github.com/shihao-jiang/meltinv.git
cd meltinv
python3.11 -m venv venv
source venv/bin/activate
pip install .
```

This also installs all dependencies, including a patched version of [meltPT](https://github.com/shihao-jiang/meltPT_fixed) that is used for the fractionation correction. If you want to modify the package, install it in editable mode instead:

```bash
pip install -e .
```

Check that the command-line tool is available:

```bash
meltinv --help
```

---

## 2. Tutorial: your first inversion

The repository ships with **`examples.xlsx`**, which contains 38 basalts from five Austral–Cook islands (Mangaia, Raivavae, Rapa, Atiu and Marotiri). We will use it throughout this tutorial.

> All outputs are written to `inversion_results/` and `inversion_figures/` inside your **current working directory**, not next to the input file. Run meltinv from the folder where you want the results to appear.

### Step 1: Prepare the input file

meltinv reads **`.xlsx`** or **`.csv`** files with one row per sample. Extra columns are allowed and are carried through to the outputs.

**Always required:**

| Column | Meaning |
|---|---|
| `location` | Group label. Samples with the same `location` are inverted together. |
| `La`, `Ce`, `Pr`, `Nd`, `Sm`, `Eu`, `Gd`, `Tb`, `Dy`, `Ho`, `Er`, `Tm`, `Yb`, `Lu` | REE concentrations (ppm) |

**Required only if you run the fractionation correction (Step 2):**

| Column | Meaning |
|---|---|
| `SiO2`, `TiO2`, `Al2O3`, `FeO`, `CaO`, `MgO`, `MnO`, `K2O`, `Na2O`, `P2O5` | Major-element oxides (wt%) of the measured (evolved) sample |

**Optional:**

| Column | Meaning |
|---|---|
| `thickness` | An independent estimate of lithospheric thickness (km), for example from seismology. It is **not** used in the inversion. Its group mean is reported as `litho_thickness_ref` and drawn as a dashed line on the figures for comparison. |

Tips for the REE columns:

- Missing values are fine. Blank cells, `0`, negative values and text (such as `bdl`) are all treated as missing.
- In group mode, an element is only used for a location if at least **2** samples in that group have a value for it. In per-sample mode the limit is 1.

### Step 2: (Optional) Fractionation correction

The inversion compares your data with **primary-magma** REE concentrations. If your samples are evolved, first correct them back to a primary composition. meltinv does this with meltPT, which adds olivine incrementally until the melt is in equilibrium with mantle olivine of forsterite content `src_Fo`.

**Command line:**

```bash
meltinv correction examples.xlsx --src_Fo 0.9 --max_olivine_addition 0.5
```

**Python:**

```python
import meltinv

corrected = meltinv.correction("examples.xlsx", src_Fo=0.9, max_olivine_addition=0.5)
```

This writes `inversion_results/examples_corrected.xlsx`. The file contains:

- all your original columns, with the **REE columns replaced** by primary-magma values (diluted by the added olivine, `REE × (1 − ol_added)`);
- `ol_added`, the mass fraction of olivine added to each sample;
- `SiO2_primary`, `TiO2_primary`, …, `P2O5_primary`, the reconstructed primary major-element composition.

> If the input is an Excel file, a `.csv` copy is also written next to it, because meltPT reads CSV.

**Skip this step** if your file already contains primary-magma REE concentrations.

### Step 3: Run the inversion

You can either invert a file that is already corrected, or correct and invert in one go.

**Option A: correct and invert in one command** (recommended for the example):

```bash
meltinv invert examples.xlsx --correction --individual_sample_figures
```

```python
import meltinv

results = meltinv.invert("examples.xlsx", correction=True, individual_sample_figures=True)
```

**Option B: invert data that is already primary**, such as the output of Step 2:

```bash
meltinv invert inversion_results/examples_corrected.xlsx --individual_sample_figures
```

```python
results = meltinv.invert("inversion_results/examples_corrected.xlsx", individual_sample_figures=True)
```

While it runs, meltinv prints each location (in alphabetical order) and its sample count, followed by every enrichment value that improves the fit. With `--correction`, meltPT first prints a long log of the olivine addition. The example takes under a minute.

```text
Austral-Cook Islands (Atiu) Sample number: 3
Scale: 0.0 Global min misfit: ...
Scale: 0.1 Global min misfit: ...
...
Figure saved to: .../inversion_figures/ree_comparison_Austral-Cook Islands (Atiu)_log10_N3.png
Austral-Cook Islands (Mangaia) Sample number: 14
...
Results saved to: .../inversion_results/examples_group_inversion_summary_log10.xlsx
```

Notice that **Marotiri is missing from the output**. It has only one sample, and group mode needs **at least 2 samples per location**. To include single samples, see [per-sample mode](#inverting-each-sample-individually).

### Step 4: Reading the results

The summary is saved to `inversion_results/examples_group_inversion_summary_log10.xlsx`. There is one row per location. For the example, the main columns look like this:

| location | count | T_mean_cali (K) | depth_mean_cali (km) | basalt_percentage | litho_thickness_ref (km) |
|---|---|---|---|---|---|
| Austral-Cook Islands (Atiu) | 3 | 70.5 | 58.7 | 4.0 | 86.7 |
| Austral-Cook Islands (Mangaia) | 14 | 91.8 | 54.2 | 1.1 | 83.4 |
| Austral-Cook Islands (Raivavae) | 14 | 92.6 | 50.7 | 0.3 | 73.7 |
| Austral-Cook Islands (Rapa) | 6 | 172.2 | 94.3 | 4.6 | 67.9 |

All output columns:

| Column | Meaning |
|---|---|
| `location`, `count` | Group name and number of samples in it |
| `T_mean` | Best-fit mantle excess temperature (K) from the peridotite melting grid |
| `T_mean_cali` | `T_mean` after the source-lithology calibration (see [§5](#5-how-the-inversion-works)). **Use this value.** |
| `depth_mean` | Best-fit lithospheric thickness (km) from the grid |
| `depth_mean_cali` | Calibrated lithospheric thickness (km). **Use this value.** |
| `pres_mean_cali`, `pres_min_cali`, `pres_max_cali` | Mean, minimum and maximum melting depth (km, about 33 km/GPa), shifted by the same calibration as the thickness |
| `F_mean`, `F_min`, `F_max` | Mean, minimum and maximum melt fraction of the best-fit model |
| `M_mean`, `M_min`, `M_max` | Mean, minimum and maximum melting rate of the best-fit model |
| `basalt_percentage` | Best-fit **relative source enrichment**. Positive values (0 to 15) are the % of basalt mixed into primitive mantle. Negative values (−100 to −5, only with `--depleted`) are the % of depleted mantle mixed into primitive mantle. |
| `misfit_value` | Misfit of the best-fit model (lower is better) |
| `global_min_value` | Lowest misfit anywhere in the full temperature–thickness–enrichment search |
| `ratio_of_local_global_misfit` | `misfit_value / global_min_value`. A value of 1 means the reported solution is the global minimum. |
| `litho_thickness_ref` | Mean of your optional `thickness` column (empty if not supplied) |

### Step 5: Reading the figures

meltinv can make two kinds of figure. Both are saved in `inversion_figures/`.

#### Individual figures (`--individual_sample_figures`)

With `--individual_sample_figures`, one PNG per location is saved to `inversion_figures/ree_comparison_{location}_log10_N{count}.png`. Each figure has four panels:

| Panel | What it shows |
|---|---|
| **Top left** | Chondrite-normalised REE patterns: samples (black), the best-fit model (red) with its uncertainty envelope (shading), and the global best fit (orange). The inverted temperature, thickness and enrichment are printed in the corner. |
| **Top right** | Misfit across excess temperature and lithospheric thickness, at the best-fit enrichment. Colours and contours show misfit divided by the minimum misfit, so 1.1 means 10% worse than the best fit. The red star marks the best fit, the orange circle the global minimum, and the dashed line `litho_thickness_ref`. |
| **Bottom left** | Misfit across excess temperature and enrichment, at the best-fit thickness |
| **Bottom right** | Misfit across thickness and enrichment, at the best-fit temperature |

The contours give a feel for how tightly each parameter is constrained. A long, narrow low-misfit valley, as in the temperature–thickness panel, shows a trade-off between the two parameters.

---

## 3. Going further

### Depleted mantle sources

By default, enrichment is searched from 0 to 15% basalt added to primitive mantle. Some locations, such as Iceland or the Galápagos, may instead sample a source **more depleted** than primitive mantle. Use `--depleted` to also search mixtures of primitive mantle with depleted mantle (−100 to −5%) for those locations:

```bash
# Locations whose name contains "Iceland" or "Galapagos"
meltinv invert input.xlsx --depleted Iceland Galapagos

# Every location
meltinv invert input.xlsx --depleted all
```

```python
meltinv.invert("input.xlsx", depleted_locations=["Iceland", "Galapagos"])
meltinv.invert("input.xlsx", depleted_locations="all")
```

Matching is by **substring and case-sensitive**: `Iceland` matches `Iceland (NNVZ)` and `Iceland (WVZ)`, but not `iceland`.

### Inverting each sample individually

Group mode fits one solution to all samples at a location. Per-sample mode inverts **each sample on its own**. This lets you see the spread within a location, and it also includes locations that have only one sample:

```bash
meltinv invert examples.xlsx --correction --per_sample
```

```python
meltinv.invert("examples.xlsx", correction=True, per_sample=True)
```

This writes two files:

| File | Contents |
|---|---|
| `inversion_results/{name}_sample_inversion_log10.xlsx` | One row per sample, in input order, with an extra `sample_number` column (1, 2, … within each location) |
| `inversion_results/{name}_sample_inversion_summary_log10.xlsx` | One row per location: a misfit-weighted average of the per-sample results (see below) |

How the location summary is built:

1. For locations with fewer than 6 samples, all samples are kept. Otherwise, the best 50% by misfit are kept, and never fewer than 3.
2. Every column is averaged with weights of `1 / misfit_value`.
3. `T_mean_cali`, `depth_mean_cali` and `basalt_percentage` also get `_err_minus` and `_err_plus` columns (weighted mean ± weighted standard deviation). For `basalt_percentage`, these error bounds use a compressed scale in which negative (depleted) values are divided by 100.

#### Group comparison figures

In per-sample mode, meltinv also makes one **group comparison figure** per location by default. It is saved as `inversion_figures/group_comparison_{location}_N{count}.png` and compares the per-sample results within that location:

| Panel | What it shows |
|---|---|
| **(a)** | Excess temperature vs lithospheric thickness, one point per sample |
| **(b)** | Relative enrichment vs lithospheric thickness |
| **(c)** | Excess temperature vs relative enrichment |
| **(d-1) to (d-3)** | Histograms of temperature, thickness and enrichment. The dashed line is the misfit-weighted mean and the grey band is ± 1 weighted standard deviation, computed the same way as the location summary above. |

In panels (a) to (c), points are coloured by misfit; darker blue is a better fit. Depleted enrichment values (negative) are divided by 100 so they fit on the same axis. Locations with fewer than 3 samples are skipped, with a message. The axes are fixed (temperature 50 to 300 K, thickness 15 to 130 km, enrichment −1.2 to 16) so that figures for different locations can be compared directly; samples outside these ranges are not shown.

To turn this figure off, use `--no_group_comparison_figure` (Python: `group_comparison_figure=False`).

With `--individual_sample_figures`, per-sample mode also makes one individual figure per sample, named `ree_comparison_{location}__sample{n}_log10_N1.png`. For large datasets this can be many files.

### Using the results in Python

`meltinv.invert` returns the results table as a pandas DataFrame, so you can analyse it directly:

```python
import meltinv

df = meltinv.invert("examples.xlsx", correction=True)

print(df[["location", "T_mean_cali", "depth_mean_cali", "basalt_percentage"]])

# Compare inverted thickness with an independent estimate
df["thickness_diff"] = df["depth_mean_cali"] - df["litho_thickness_ref"]
```

In per-sample mode, the returned DataFrame holds the per-sample results. The location summary is only written to file.

---

## 4. Reference

### Python API

```python
meltinv.correction(file_name, src_Fo=0.9, max_olivine_addition=0.5)
```

| Argument | Default | Description |
|---|---|---|
| `file_name` | — | Path to `.xlsx` or `.csv` input |
| `src_Fo` | `0.9` | Forsterite content of mantle-source olivine |
| `max_olivine_addition` | `0.5` | Maximum mass fraction of olivine that may be added |

Returns the corrected DataFrame and writes `inversion_results/{name}_corrected.xlsx`.

```python
meltinv.invert(file_name, depleted_locations=None, correction=False, src_Fo=0.9,
               max_olivine_addition=0.5, individual_sample_figures=False,
               group_comparison_figure=True, per_sample=False)
```

| Argument | Default | Description |
|---|---|---|
| `file_name` | — | Path to `.xlsx` or `.csv` input |
| `depleted_locations` | `None` | List of location substrings, or `"all"`, for which depleted sources are also searched |
| `correction` | `False` | Run the fractionation correction first |
| `src_Fo` | `0.9` | Used only when `correction=True` |
| `max_olivine_addition` | `0.5` | Used only when `correction=True` |
| `individual_sample_figures` | `False` | Save a 4-panel REE and misfit figure for each group (group mode) or each sample (per-sample mode). The old name `make_figures` still works. |
| `group_comparison_figure` | `True` | In per-sample mode, save a group comparison figure for each location with ≥ 3 samples. Ignored in group mode. |
| `per_sample` | `False` | Invert each sample individually instead of each location |

Returns the results as a pandas DataFrame.

### Command line

```text
meltinv correction FILE [--src_Fo 0.9] [--max_olivine_addition 0.5]

meltinv invert FILE [--depleted LOC [LOC ...] | --depleted all]
                    [--correction] [--src_Fo 0.9] [--max_olivine_addition 0.5]
                    [--individual_sample_figures] [--no_group_comparison_figure]
                    [--per_sample]
```

The command-line flags map onto the Python arguments above. `--no_group_comparison_figure` is the same as `group_comparison_figure=False`, and the old `--make_figures` flag still works as `--individual_sample_figures`.

### Output files

All paths are relative to the current working directory. `{name}` is the input file name without its extension; a trailing `_corrected` is dropped.

| File | Produced by |
|---|---|
| `inversion_results/{name}_corrected.xlsx` | `correction`, or `invert --correction` |
| `inversion_results/{name}_group_inversion_summary_log10.xlsx` | `invert` (group mode) |
| `inversion_results/{name}_sample_inversion_log10.xlsx` | `invert --per_sample` |
| `inversion_results/{name}_sample_inversion_summary_log10.xlsx` | `invert --per_sample` |
| `inversion_figures/ree_comparison_{location}_log10_N{count}.png` | `invert --individual_sample_figures` |
| `inversion_figures/group_comparison_{location}_N{count}.png` | `invert --per_sample` (unless `--no_group_comparison_figure`) |

Existing files with the same name are overwritten.

---

## 5. How the inversion works

1. **Forward-model grid.** The package ships with precomputed melting models (`meltinv/data/*.npy`). They give the REE concentrations of accumulated melts, with min/max envelopes, melting depths, melt fractions (`F`) and melting rates (`M`). The grid spans mantle excess temperature (65 to ~350 K in 5 K steps) and lithospheric thickness (15 to ~120 km in 2 km steps).
2. **Source enrichment.** For each trial enrichment value, the model REE concentrations are rescaled for a source made by mixing:
   - primitive mantle (McDonough & Sun, 1995) with 0 to 15% basalt (Sun & McDonough, 1989), in steps of 0.1%; or
   - with `--depleted`, primitive mantle with 5 to 100% depleted mantle (Salters & Stracke, 2004).
3. **Misfit.** For each element, meltinv computes the log10 misfit between every sample and the model, takes the root-sum-square, and divides by the number of samples. It then averages over all usable elements.
4. **Best fit.** The temperature, thickness and enrichment combination with the lowest misfit is reported.
5. **Calibration.** The grid assumes a peridotite source. The best-fit temperature and thickness are then corrected for the effect of a pyroxenite (G2) component on the solidus. This uses flux-weighted mixing (Phipps Morgan, 2001; Sobolev et al., 2007), with the pyroxenite fraction taken from the inverted enrichment. The results are the `_cali` columns.

---

## 6. Troubleshooting

| Problem | Cause and fix |
|---|---|
| `The following required columns are missing: ...` | Column names are case-sensitive. Check `location` and the 14 REE names. |
| `Unsupported file format` | Only `.xlsx`, `.xls` and `.csv` are accepted. |
| A location is missing from the results | It has only 1 sample (group mode needs ≥ 2), or no element has values in enough samples. Use `--per_sample` to include single samples. |
| `KeyError` on `SiO2` (or another oxide) during correction | The correction needs all 10 major-oxide columns listed in [Step 1](#step-1-prepare-the-input-file). |
| Results appear in an unexpected folder | Outputs go to the **current working directory**. `cd` to where you want them before running. |

---

## 7. Citation

If you use **meltinv** in your research, please cite:

Jiang et al. (2026), xxx, *Journal Name*.

## License

MIT. See [LICENSE](LICENSE).
