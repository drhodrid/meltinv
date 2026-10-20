import pandas as pd
import numpy as np
import warnings
from pathlib import Path

def get_and_check_file(file_name):
    '''
    Read melt composition data from an Excel or CSV file.

    The input file must contain the required rare earth element (REE) concentrations in the primary magma
    (i.e., fractionation corrected)along with a sample location column. This function validates the file format
    and checks that all required columns are present before returning the data as a pandas DataFrame.
    '''
    file_path = Path(file_name).expanduser().resolve()

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_name}.")
    suffix = file_path.suffix.lower()

    if suffix in [".xlsx", ".xls"]:
        df = pd.read_excel(file_path, engine="openpyxl")
        new_path = file_path.with_suffix(".csv")
        df.to_csv(new_path, index=False)
        print(f"Converted Excel file to CSV: {new_path}")
        file_path = new_path
    elif suffix in [".csv"]:
        df = pd.read_csv(file_path)
    else:
        raise ValueError("Unsupported file format. Please provide a .xlsx or .csv file.")

    required_columns = ['location',
                        'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd',
                        'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']

    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        missing_str = ", ".join(missing)
        raise ValueError(f"The following required columns are missing: {missing_str}")
    else:
        print("Successfully retrieved the sample composition file for fractionation correction.")

    return df, file_path

def save_fraction(file_name, sample_df, primary_composition_df):
    results_dir = Path("inversion_results")
    results_dir.mkdir(exist_ok=True)

    base_name = Path(file_name).stem
    output_file = results_dir / f"{base_name}_corrected.xlsx"

    # A correction counts as successful only if meltPT returned a valid olivine addition
    # and a complete primary major-element composition; otherwise the REE are set to missing
    sample_df = sample_df.copy()
    primary_composition_df = primary_composition_df.reindex(sample_df.index)
    olv_added = pd.to_numeric(primary_composition_df['ol_added'], errors='coerce')
    primary_cols = [f'{col}_primary_wt' for col in
                    ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'CaO', 'MgO', 'MnO', 'K2O', 'Na2O', 'P2O5']]
    valid = (np.isfinite(olv_added) & olv_added.between(0, 1, inclusive='left')
             & np.isfinite(primary_composition_df[primary_cols]).all(axis=1))
    if not valid.all():
        warnings.warn(f"Fractionation correction failed for {(~valid).sum()} sample(s); "
                      "their corrected REE values are missing", RuntimeWarning)
    olv_added = olv_added.where(valid)
    sample_df["correction_success"] = valid
    columns_to_multiply = ['La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    sample_df["ol_added"] = olv_added
    sample_df[columns_to_multiply] = (sample_df[columns_to_multiply].apply(pd.to_numeric, errors="coerce")
                                      .mul(1 - olv_added, axis=0))

    major_columns = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'CaO', 'MgO', 'MnO', 'K2O', 'Na2O', 'P2O5']

    for col in major_columns:
        insert_at = sample_df.columns.get_loc(col) + 1
        sample_df.insert(insert_at, f'{col}_primary',
                         primary_composition_df[f'{col}_primary_wt'])

    sample_df.to_excel(output_file, index=False)
    print(f"Results saved to: {output_file.resolve()}")

    return sample_df

def fractionation_correction(file_name, src_Fo=0.9, max_olivine_addition=0.4):
    # Import meltPT only when a correction is requested
    from meltPT import Suite, BacktrackOlivineFractionation

    sample_df, file_path = get_and_check_file(file_name)
    s = Suite(file_path, src_Fo=src_Fo)
    b = BacktrackOlivineFractionation(verbose=True, max_olivine_addition=max_olivine_addition)
    s.backtrack_compositions(backtracker=b)
    primary_composition_df = s.primary
    output_df = save_fraction(file_name, sample_df, primary_composition_df)

    return output_df




