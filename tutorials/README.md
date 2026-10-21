# meltinv tutorials

Step-by-step Jupyter notebooks that take you through a meltinv analysis. You can read them online, run them in your browser without installing anything, or run them on your own computer.

[![Launch on Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/shihao-jiang/meltinv/main?labpath=tutorials%2F01_quick_start.ipynb)

| Tutorial | What you will do | Data | Run time |
|---|---|---|---|
| [01 – Quick start](01_quick_start.ipynb) | Correct, invert and interpret the bundled example: input checks, olivine correction, results tables, figures and a comparison with seismic thickness | `examples.xlsx` (38 samples) | ~1 min |
| 02 – Modern ocean-island basalts *(in preparation)* | The full dataset used in the paper, ending with the comparison figure | modern OIB (1,652 samples) | — |

## Reading the tutorials online

The notebooks are saved with all their outputs and figures, so you can read them without running anything:

- on GitHub, by clicking a notebook in the table above;
- on [nbviewer](https://nbviewer.org/github/shihao-jiang/meltinv/blob/main/tutorials/01_quick_start.ipynb), which renders notebooks a little more faithfully.

## Running the tutorials in your browser (Binder)

Click the **launch binder** button above. [Binder](https://mybinder.org) builds a temporary JupyterLab with meltinv installed and opens Tutorial 1. The first launch can take a few minutes. Binder sessions are temporary, so download anything you want to keep before closing the tab.

## Running the tutorials on your own computer

Install meltinv with the notebook tools, from the repository folder:

```bash
pip install -e ".[tutorials]"
```

Then start JupyterLab from the `tutorials` folder, using the same Python that meltinv is installed in:

```bash
cd tutorials
python -m jupyterlab
```

Open a notebook and run the cells from top to bottom (**Run → Run All Cells**).

If the first cell fails with `ModuleNotFoundError: No module named 'meltinv'`, the notebook is running in a different Python from the one meltinv is installed in. Close JupyterLab and start it again with `python -m jupyterlab` from the environment where you installed meltinv. On macOS, plain `jupyter lab` can pick the wrong Python.

Each tutorial writes its outputs to its own folder, `tutorials/output/<tutorial name>/`, so nothing is written into the rest of the repository. Delete that folder at any time to start again.
