## Code to add plots to axes for correction factor.
## This is just to save the crazy amount of code already used to generate Figure 3.

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D

import statsmodels.formula.api as smf


from paiskintonetools.stats import loess_bootstrap  # type: ignore


def p_formater(p):
    if p < 0.001:
        return "$p < 0.001$"
    else:
        return f"$p = {p:.3f}$"


def agg_function(x):
    return (
        np.mean(np.stack(x), axis=0) if np.issubdtype(x.dtype, np.number) else x.iloc[0]
    )


def mua_melanin(x, mvf):
    return 519 * (x / 500) ** (-3.5) * mvf / 10  # /mm


def make_figure_6_plot(ax, ax_cf, axso2, colours):
    df = pd.read_csv("../Fluence Correction/tables/cali_curve.csv", index_col=0)
    df.set_index(["MVF", "WL"], inplace=True, drop=False)

    df["Compare Fluence"] = np.exp(-mua_melanin(df["WL"], df["MVF"]) * 0.06)

    for n, g in df.groupby(level=0):
        # print(df.loc[(0.02, ), "Fluence"])
        df.loc[(n,), "Normalised"] = (
            df.loc[(n,), "Fluence"].values / df.loc[(0.02,), "Fluence"].values  # type: ignore
        )
        df.loc[(n,), "Normalised Baseline"] = (
            df.loc[(n,), "Compare Fluence"].values  # type: ignore
            / df.loc[(0.02,), "Compare Fluence"].values  # type: ignore
        )
    all_mvfs = pd.unique(df["MVF"])

    ITA = np.array([67, 56, 41, 20, -6.9, -47])
    MVF = np.array([2, 3.6, 6.6, 12.1, 22, 40])
    ax.scatter(ITA, MVF, s=5, c="k")
    ax.set_xlabel("ITA (degrees)")
    ax.set_ylabel("Melanosome volume\nfraction (%)")
    ols_fit = smf.ols(
        "MVF ~ ITA + np.power(ITA, 2)", data=pd.DataFrame({"ITA": ITA, "MVF": MVF})
    )

    f = ols_fit.fit()
    ita_eval = pd.DataFrame({"ITA": np.linspace(np.min(ITA), np.max(ITA), 100)})
    p = f.get_prediction(ita_eval)
    ita_eval = pd.concat([ita_eval, p.summary_frame()], axis=1)

    constant = f.params["Intercept"]
    ita_coeff = f.params["ITA"]
    ita2_coeff = f.params["np.power(ITA, 2)"]
    label = f"{constant:.0f} $-$ {-ita_coeff:.2f} $\\times$ ITA\n + {ita2_coeff:.4f} $\\times$ ITA$^2$"

    ax.plot(ita_eval["ITA"], ita_eval["mean"], label=label)
    ax.invert_xaxis()

    ax.legend(fontsize="small", loc="lower right")
    ax_cf.set_ylabel("log(Correction factor)")
    ax_cf.set_xlabel("Wavelength (nm)")

    mvf_max = np.max(df["MVF"] * 100)
    mvf_min = np.min(df["MVF"] * 100)

    norm = mpl.colors.LogNorm(vmin=1, vmax=mvf_max)  # type: ignore
    scm = plt.cm.ScalarMappable(cmap="YlGnBu", norm=norm)

    i = 0
    for n, g in df.groupby(level=0):
        plot_line = ax_cf.plot(
            np.linspace(700, 900, 5),
            g["Normalised"],
            color=scm.to_rgba(n * 100),  # type: ignore
        )
        ax_cf.plot(
            np.linspace(700, 900, 5),
            g["Normalised Baseline"],
            color=scm.to_rgba(n * 100),  # type: ignore
            linestyle="--",
        )
        i += 1

    l0 = Line2D([0], [0], linewidth=plot_line[0].get_linewidth(), color="k")  # type: ignore
    l1 = Line2D(
        [0],
        [0],
        linewidth=plot_line[0].get_linewidth(),  # type: ignore
        color="k",
        linestyle="--",  # type: ignore
    )
    ax_cf.legend([l0, l1], ["Monte-Carlo", "Beer-Lambert"], fontsize="small")

    cbar = plt.colorbar(scm, ax=ax_cf, aspect=15)
    cbar.set_label(
        "Melanosome volume\nfraction (%)", fontsize="small", fontweight="normal"
    )
    cbar.ax.set_yticks(
        100 * all_mvfs,  # type: ignore
        [f"{100 * x:.1f}" for i, x in enumerate(all_mvfs)],
    )
    cbar.ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())  # type: ignore
    cbar.ax.set_ylim([mvf_min, mvf_max])  # type: ignore
    # f = mpl.ticker.ScalarFormatter()
    # f.set_scientific(True)
    # cbar.ax.yaxis.set_major_formatter(f)

    y_plot = "so2_mean"

    df_data = pd.read_parquet("intermediate data/forearm_artery_results_all.parquet")
    df_data = (
        df_data.groupby(["SkinID", "Region", "RepNumber"])
        .agg(agg_function)
        .groupby(level=[0, 1])
        .agg(agg_function)
        .reset_index()
    )

    df_data_bicep = pd.read_parquet(
        "intermediate data/bicep_muscle_results_all.parquet"
    )
    df_data_bicep = (
        df_data_bicep.groupby(["SkinID", "Region", "RepNumber"])
        .agg(agg_function)
        .groupby(level=[0, 1])
        .agg(agg_function)
        .reset_index()
    )

    df_ita = (
        pd.read_parquet("../SummaryTables/ita_etc.parquet")
        .loc[(slice(None, None), "Forearm"), slice(None, None)]
        .reset_index()
    )
    df_fp = pd.read_excel("../SummaryTables/SummaryDetails.xlsx")
    df_ita = df_ita.merge(df_fp, on="SkinID")

    df_data = df_data.merge(df_ita, on="SkinID")

    df_ita_bicep = (
        pd.read_parquet("../SummaryTables/ita_etc.parquet")
        .loc[(slice(None, None), "Bicep"), slice(None, None)]
        .reset_index()
    )
    df_ita_bicep = df_ita_bicep.merge(df_fp, on="SkinID")
    df_data_bicep = df_data_bicep.merge(df_ita_bicep, on="SkinID")

    for region, col in zip(["forearm", "bicep"], colours):
        df_toplot = df_data_bicep if region == "bicep" else df_data

        # ax.scatter(df["ITA"], df[y_plot], s=5)
        X, Y, Ymax, Ymin = loess_bootstrap(
            df_toplot["ITA"], df_toplot[y_plot], frac=0.6666, it=1, seed=1
        )
        plot_line = axso2.plot(X, Y * 100, label="Uncorrected", c=col)
        axso2.axhline(Y[-1] * 100, c="gray", linewidth=1)
        # axso2.fill_between(X, Ymin, Ymax, alpha=0.3, facecolor=l[0].get_color())

        # ax.scatter(df["ITA"], df["corrected_" + y_plot], s=5)
        X, Y, Ymax, Ymin = loess_bootstrap(
            df_toplot["ITA"],
            df_toplot["corrected_" + y_plot],
            frac=0.66666,
            it=1,
            seed=1,
        )
        plot_line = axso2.plot(X, Y * 100, label="Corrected", c=col, linestyle="--")
        axso2.fill_between(
            X, Ymin * 100, Ymax * 100, alpha=0.3, facecolor=plot_line[0].get_color()
        )
    axso2.invert_xaxis()

    handles, labels = axso2.get_legend_handles_labels()
    handles = [
        Line2D(
            [0],
            [0],
            linewidth=handle.get_linewidth(),
            color="k",
            linestyle=handle.get_linestyle(),
        )
        for handle in handles[:2]
    ]
    labels = labels[:2]
    l0 = Line2D([0], [0], linewidth=plot_line[0].get_linewidth(), color=colours[0])  # type: ignore
    l1 = Line2D([0], [0], linewidth=plot_line[0].get_linewidth(), color=colours[1])  # type: ignore
    handles += [l0, l1]
    labels += ["Radial artery", "Bicep muscle"]
    axso2.legend(handles, labels, fontsize="small")
    axso2.set_xlabel("ITA (degrees)")
    axso2.set_ylabel(r"sO$_\mathbf{2}^\mathbf{EST}$ (%)")
