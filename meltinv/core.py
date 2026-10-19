from . import inversion
from .inversion import invert_melt_condition
from .frac_correction import fractionation_correction

def invert(file_name, depleted_locations=None, correction=False, src_Fo=0.9,
           max_olivine_addition=0.5, individual_sample_figures=False,
           group_comparison_figure=True, per_sample=False, make_figures=None):
    """
    High-level API for melt inversion.

    make_figures is kept for backward compatibility; it sets individual_sample_figures.
    """
    if make_figures is not None:
        individual_sample_figures = make_figures

    summary_df = invert_melt_condition(file_name, depleted_locations, correction,
                          src_Fo, max_olivine_addition, individual_sample_figures,
                          group_comparison_figure, per_sample)

    return summary_df

def correction(file_name, src_Fo=0.9, max_olivine_addition=0.5):
    return fractionation_correction(file_name, src_Fo, max_olivine_addition)



