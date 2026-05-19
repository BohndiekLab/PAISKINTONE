from patato.recon.fourier_transform_rec import FFTReconstruction
import patato as pat
from scipy.signal import butter, filtfilt
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.sparse import load_npz
from scipy.interpolate import interp1d

from tqdm.auto import tqdm

from pathlib import Path


def ft_to_image(fd):
    return np.fft.fftshift(np.fft.fft2(np.fft.fftshift(fd)))


def reconstruct_plane_wave_us(ft, dx_move, dx, return_ft=False, scanner_move=0):
    ft_data = np.fft.fftshift(ft.T)
    kz = kx = np.fft.fftfreq(ft_data.shape[0])
    kz = kz[:, None]
    kx = kx[None, :]
    kz_primes = np.sqrt(kx**2 + kz**2) + kz

    ft_data = ft_data * np.exp(1j * 2 * np.pi * kz * dx_move / dx)

    ongrid = np.ones((kx.shape[1] * 3, kx.shape[1]), dtype=np.complex128)

    kz_interp = np.linspace(kx.min(), kx.max(), kx.shape[1] * 3)

    for i in range(kx.shape[1]):
        kz_prime = kz_primes[:, i]
        nonans = ~np.isnan(kz_prime)
        ongrid[:, i] = interp1d(
            kz_prime[nonans],
            np.real(ft_data[:, i])[nonans],
            fill_value=0,
            bounds_error=False,
        )(kz_interp)
        ongrid[:, i] += 1j * interp1d(
            kz_prime[nonans],
            np.imag(ft_data[:, i])[nonans],
            fill_value=0,
            bounds_error=False,
        )(kz_interp)
    ongrid *= np.exp(-1j * 2 * np.pi * kz_interp * dx_move / dx)[:, None]
    ongrid *= np.exp(1j * 2 * np.pi * kz_interp * scanner_move / dx)[:, None]
    ongrid = np.fft.fftshift(ongrid, 1)
    if return_ft:
        return ongrid, ft_to_image(ongrid)
    return ft_to_image(ongrid)


if __name__ == "__main__":
    matplotlib.use("Agg")
    dx_move = 0.00474
    dx = 0.04 / 399

    import json
    from pathlib import Path

    settings = json.load(open("../data_paths.json", "r"))
    data_path = (Path.cwd() / "../" / settings["main_data_path"]).resolve()
    us_path = (Path.cwd() / "../" / settings["melanin_ultrasound_path"]).resolve()

    data_folder = data_path
    output_folder = us_path.parent

    # Set up reconstruction class
    fft_rec = FFTReconstruction(field_of_view=(0.04, 0, 0.04), n_pixels=(400, 1, 400))
    fft_rec.hankels = None

    # Get model matrix
    m = load_npz(settings["model_path"])

    files_to_process = list(data_folder.glob("**/*.hdf5"))

    for pa_file in tqdm(files_to_process, desc="Scan loop"):
        pa = pat.PAData.from_hdf5(pa_file)
        name = pa.get_scan_name()
        if "3" in name or "2" in name:
            continue
        scan_id = pa_file.stem
        skin_id = pa_file.parent.stem

        output_folder_id = output_folder / skin_id
        output_folder_id.mkdir(parents=True, exist_ok=True)

        if ("skin_", "0") not in pa.get_rois():
            print(f"Skipping {skin_id}, {scan_id}, {name}.")
            continue

        rec = pa.get_scan_reconstructions()["Model Based", "0"]
        rec_raw_data = rec.raw_data

        rec_data = np.squeeze(rec_raw_data)  # type: ignore
        mask = np.squeeze(pa.get_rois()["skin_", "0"].to_mask_slice(rec)[0])

        timeseries = np.squeeze(pa.get_time_series().raw_data)  # type: ignore

        us_images = []
        us_images_filtered = []
        for i in tqdm(range(rec_data.shape[0]), leave=False, desc="Wavelength loop"):
            rec_data[i][(~mask) | (rec_data[i] < 0)] = 0  # type: ignore
            ts_melanin = (m @ rec_data[i].ravel()).reshape((256, 2030))

            b, a = butter(5, [5e3, 7e6], btype="bandpass", fs=4e7)  # type: ignore
            ts_for_us_rec = filtfilt(
                b, a, timeseries[i], axis=-1, padtype="even", padlen=100, method="pad"
            )

            ts_for_us_rec /= pa.get_overall_correction_factor()[0, i]
            ts_for_us_rec -= ts_melanin

            _, ft, _ = fft_rec._reconstruct(
                ts_for_us_rec,
                pa.get_sampling_frequency(),
                pa.get_scan_geometry(),
                (400, 1, 400),
                (0.04, 0, 0.04),
                1520,
                return_ft=True,
                debug=False,
            )
            us_image = reconstruct_plane_wave_us(ft, dx_move, dx)
            us_images.append(us_image)
            us_image = np.abs(us_image)
            us_image = us_image[1600 - 40 : 2000 + 40, 360:840]

            # Weight and apply log compression
            weighting = np.exp(np.arange(480) / 75)
            to_display = np.abs(us_image * weighting[:, None])
            to_display = np.log10(1 + 3 * to_display) / np.log10(4)
            us_images_filtered.append(to_display)
            fig, (ax1, ax2) = plt.subplots(
                1, 2, figsize=(10, 4), sharex=True, sharey=True
            )
            im = ax1.imshow(
                to_display, extent=(-0.02 * 1.1, 0.02 * 1.1) * 2, cmap="magma"
            )
            plt.colorbar(im, ax=ax1)
            im = ax2.imshow(
                pa.get_ultrasound().raw_data[0, 0][:, 0],  # type: ignore
                extent=(-0.02, 0.02) * 2,
                origin="lower",
                cmap="bone",
            )
            plt.colorbar(im, ax=ax2)
            ax1.axis("off")
            ax2.axis("off")
            ax1.set_title("Optical ultrasound")
            ax2.set_title("Ultrasound")
            fig.suptitle(f"{skin_id} {scan_id} {name}")
            fig.tight_layout()
            fig.savefig(
                output_folder_id / (str.zfill(str(i), 2) + "_" + scan_id + "_plot"),
                dpi=300,
            )
            plt.close()
        output_images = np.stack(us_images_filtered)
        output_us = np.stack(us_images)
        np.save(output_folder_id / (scan_id + "_images"), output_images)
        np.save(output_folder_id / (scan_id + "_raw_us"), output_us)
