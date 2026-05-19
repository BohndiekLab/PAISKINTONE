import numpy as np
import patato as pat
from patato.io.msot_data import HDF5Tags
from patato.io.attribute_tags import ReconAttributeTags

from pathlib import Path

import cupy as cp  # type: ignore

from scipy.signal import butter, filtfilt
from cupyx.scipy.sparse import csr_matrix  # type: ignore

from pylops.optimization.basic import lsqr
from pylops import MatrixMult

import json

settings = json.load(open(Path(__file__).parent.parent / "data_paths.json", "r"))

forward_matrix = np.load(settings["model_path"])
forward_matrix = csr_matrix(
    (
        cp.array(forward_matrix["data"]),
        cp.array(forward_matrix["indices"]),
        cp.array(forward_matrix["indptr"]),
    ),
    shape=tuple(forward_matrix["shape"]),
)
forward_matrix = MatrixMult(forward_matrix.astype(float))


root = Path("../Data")
images = root.glob("**/*.hdf5")


def solve(ts):
    b, a = butter(5, [5e3, 7e6], btype="bandpass", fs=4e7)  # type: ignore
    fp = cp.array(filtfilt(b, a, ts, axis=-1, padtype="even", padlen=100, method="pad"))
    return lsqr(forward_matrix, fp.flatten(), damp=3e2)[0].get()


def reconstruct(ts):
    rd = ts.raw_data.astype(float)
    orig_shape = rd.shape
    rd = rd.reshape((-1,) + orig_shape[-2:])
    results = []
    for i in range(rd.shape[0]):
        results.append(solve(rd[i]).reshape((400, 400)))
    return np.stack(results).reshape(orig_shape[:2] + (400, 1, 400))


def nonneg_mean(x, axis=None):
    y = np.copy(x)
    y[y < 0] = np.nan
    return np.nanmean(y, axis=axis)


def nonneg_median(x, axis=None):
    y = np.copy(x)
    y[y < 0] = np.nan
    return np.nanmedian(y, axis=axis)


def nonneg_max(x, axis=None):
    y = np.copy(x)
    y[y < 0] = np.nan
    return np.nanmax(y, axis=axis)


def pat_reconstruct(ts):
    rec_raw = reconstruct(ts)
    rec = pat.Reconstruction(
        rec_raw,
        ts.ax_1_labels,
        hdf5_sub_name="Model Based",
        field_of_view=(0.04, 0.0, 0.04),
    )
    rec.attributes[HDF5Tags.SPEED_OF_SOUND] = 1520
    rec.attributes[ReconAttributeTags.RECONSTRUCTION_ALGORITHM] = "Model Based"
    rec.attributes[ReconAttributeTags.X_NUMBER_OF_PIXELS] = 400
    rec.attributes[ReconAttributeTags.Y_NUMBER_OF_PIXELS] = 1
    rec.attributes[ReconAttributeTags.Z_NUMBER_OF_PIXELS] = 400
    rec.attributes[ReconAttributeTags.X_FIELD_OF_VIEW] = 0.04
    rec.attributes[ReconAttributeTags.Y_FIELD_OF_VIEW] = 0.0
    rec.attributes[ReconAttributeTags.Z_FIELD_OF_VIEW] = 0.04
    rec.attributes[HDF5Tags.WAVELENGTH] = ts.ax_1_labels
    rec.attributes["Damp"] = 3e2
    rec.attributes["Bandpass"] = [5e3, 7e6]
    rec.attributes["BandpassDetails"] = (
        "Scipy Butter Bandpass, 5kHz to 7MHz. Padtype even, padlen 100."
    )
    return rec


for e in images:
    print(e)
    p = pat.PAData.from_hdf5(e, "r+")  # type: ignore
    if ("Model Based", "0") in p.get_scan_reconstructions():  # type: ignore
        p.scan_reader.file.close()  # type: ignore
        p.scan_writer.file.close()  # type: ignore
        print("Skipping", e, "as reconstructed previously.")
        continue
    ts = p.get_time_series()
    rec = pat_reconstruct(ts)
    rec.raw_data /= p.get_overall_correction_factor()[:, :, None, None, None]  # type: ignore
    rec.save(p.scan_writer)
