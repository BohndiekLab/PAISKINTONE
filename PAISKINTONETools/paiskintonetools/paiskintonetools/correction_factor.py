import pandas as pd
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def get_correction_factor_interpolator(cal_curve_file):
    # Load in my fluence calibration curve
    df = pd.read_csv(cal_curve_file, index_col=0)
    df.set_index(["MVF", "WL"], inplace=True, drop=False)
    for n, g in df.groupby(level=0):
        df.loc[(n,), "Normalised"] = (
            df.loc[(n,), "Fluence"].values / df.loc[(0.02,), "Fluence"].values  # type: ignore
        )

    wavelengths = df.index.get_level_values(1)[:5].values
    mvf = df.index.get_level_values(0)[::5].values

    norm_grid = -np.log(df["Normalised"].values.reshape((-1, 5)))  # type: ignore

    # Extrapolate values outside the range we simulated. This should be fine, I think, it's a very smooth function and we're not too far outside of the range.
    correction_factor_spline = RegularGridInterpolator(
        (wavelengths, mvf),
        norm_grid.T,
        method="cubic",
        bounds_error=False,
        fill_value=None,  # type: ignore
    )
    return correction_factor_spline
