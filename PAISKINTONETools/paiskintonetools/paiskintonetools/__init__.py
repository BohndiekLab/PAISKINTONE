def setup_matplotlib(dpi=300):
    import matplotlib.pyplot as plt

    plt.rcParams["figure.dpi"] = dpi
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 8
    plt.rcParams["axes.linewidth"] = 1
    # plt.rcParams["axes.prop_cycle"] =
    plt.rcParams["legend.frameon"] = False
    plt.rcParams["lines.linewidth"] = 1
    plt.rcParams["savefig.bbox"] = None
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["xtick.major.width"] = 1
    plt.rcParams["xtick.minor.width"] = 1
    plt.rcParams["xtick.major.size"] = 3
    plt.rcParams["xtick.minor.size"] = 1.5
    plt.rcParams["xtick.minor.visible"] = True
    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["ytick.major.width"] = 1
    plt.rcParams["ytick.minor.width"] = 1
    plt.rcParams["ytick.major.size"] = 3
    plt.rcParams["ytick.minor.size"] = 1.5
    plt.rcParams["ytick.minor.visible"] = True
    plt.rcParams["ytick.right"] = True
    plt.rcParams["lines.markersize"] = 5
    plt.rcParams["axes.titleweight"] = "bold"
    plt.rcParams["axes.titlesize"] = "large"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.rcParams["figure.constrained_layout.use"] = True
    plt.rcParams["mathtext.fontset"] = "stixsans"
    plt.rcParams["svg.fonttype"] = "none"


def file_sort_key(x):
    return int(x.stem.split("_")[1])


def get_example_scan_of_type(folder, desired_type, exclude="parallel"):
    import patato as pat
    from patato.io.hdf.hdf5_interface import HDF5Reader
    from pathlib import Path
    import h5py

    for f in sorted(list(Path(folder).glob("*.hdf5")), key=file_sort_key):
        with h5py.File(f) as h5f:
            pa = pat.PAData(HDF5Reader(h5f))
            if desired_type in pa.get_scan_name() and exclude not in pa.get_scan_name():
                return f
