#!/usr/bin/env python
import matplotlib.pyplot as plt

params = {
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 14,
    "font.family": "serif",
}
plt.rcParams.update(params)
import argparse
import h5py
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("file_inv")
    parser.add_argument("--plot_model", action="store_true")
    parser.add_argument("-d", "--model_data")
    parser.add_argument("--plot_disp", action="store_true")
    parser.add_argument("--plot_fit", action="store_true")
    parser.add_argument("--full_disp", action="store_true")
    args = parser.parse_args()
    file_inv = args.file_inv
    show_model = args.plot_model
    show_disp = args.plot_disp
    show_fitness = args.plot_fit
    show_full_disp = args.full_disp
    file_model_data = args.model_data

    fh5 = h5py.File(file_inv, "r")
    z_sample = fh5["z_sample"][()]
    vs_sample = fh5["vs_sample"][()]
    vs_hist = fh5["vs_hist2d"][()]
    vs_mean = fh5["vs_mean"][()]
    vs_median = fh5["vs_median"][()]
    vs_mode = fh5["vs_mode"][()]
    vs_cred10 = fh5["vs_cred10"][()]
    vs_cred90 = fh5["vs_cred90"][()]
    mode_used = fh5["mode_used"][()]
    num_valid = fh5["num_valid"][()]
    fitness = fh5["fitness"][()]
    data = fh5["data"][()]
    disp_syn = [fh5[f"disp/{i}"][()] for i in range(num_valid)]
    fh5.close()

    if show_model:
        plot_model(
            z_sample,
            vs_sample,
            vs_hist,
            vs_mean,
            vs_median,
            vs_mode,
            vs_cred10,
            vs_cred90,
            file_model_data,
        )
    if show_disp:
        plot_disp(data, disp_syn, mode_used, show_full_disp)
    if show_fitness:
        plot_fitness(fitness)

    plt.show()


def plot_model(
    z_sample,
    vs_sample,
    vs_hist,
    vs_mean,
    vs_median,
    vs_mode,
    vs_cred10,
    vs_cred90,
    file_model_data,
):
    vs_hist = np.ma.masked_array(vs_hist, mask=vs_hist <= 0)

    fig, ax = plt.subplots(layout="constrained")
    ax.pcolormesh(
        vs_sample,
        z_sample,
        vs_hist,
        # cmap="Blues_r",
        cmap="Wistia",
        alpha=0.8,
    )

    if file_model_data:
        model_data = np.loadtxt(file_model_data)
        vs = model_data[:, 3]
        z = model_data[:, 1]
        if z[-1] < z_sample[-1]:
            z = np.append(z, z_sample[-1])
            vs = np.append(vs, vs[-1])
        ax.step(vs, z, "-", c="tab:red", alpha=0.7, lw=2, label="Target")

    ax.plot(
        vs_mean,
        z_sample,
        "-",
        c="k",
        alpha=0.7,
        lw=2,
        label="Mean",
    )
    # ax.plot(vs_mode, z_samples, "--", c="tab:red", alpha=0.7, lw=2, label="Mode")
    ax.plot(
        vs_cred10,
        z_sample,
        "k--",
        dashes=(5, 5),
        lw=1,
        alpha=0.6,
        label="10/90 percentile",
    )
    ax.plot(vs_cred90, z_sample, "k--", lw=1, alpha=0.6, dashes=(5, 5))

    ax.set_ylim([0, z_sample[-1]])
    ax.invert_yaxis()
    ax.set_xlabel("Vs (km/s)")
    ax.set_ylabel("Depth (km)")
    # ax.legend(loc="upper right")
    ax.legend()
    ax.grid(linestyle=":")


def plot_disp(data, disp_syn, mode_used, show_full_disp):
    if show_full_disp:
        mode_show = list(set(data[:, 2].astype(int)))
    else:
        mode_show = mode_used

    fig, ax = plt.subplots(layout="constrained")
    for disp in disp_syn:
        p1 = plot_1disp(ax, disp, mode_show, "", "-", "k", 0.6)
    p2 = plot_1disp(ax, data, mode_show, ".", "", "r", 0.8)
    ax.legend([p1, p2], ["inv", "data"], loc="upper right")

    modes = data[:, 2].astype(int)
    f_show = np.hstack([data[modes == m][:, 0] for m in mode_show])
    ax.set_xlim([np.amin(f_show), np.amax(f_show)])

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Phase velocity (km/s)")


def plot_1disp(ax, disp, mode_show, marker, linestyle, color, alpha):
    modes = disp[:, 2].astype(int)
    for m in mode_show:
        d = disp[modes == m]
        (p1,) = ax.plot(
            d[:, 0],
            d[:, 1],
            marker=marker,
            linestyle=linestyle,
            color=color,
            alpha=alpha,
        )
    return p1


def plot_fitness(fitness):
    idx = np.argsort(fitness)[::-1]
    fitness = fitness[idx]
    fig, ax = plt.subplots(layout="constrained")
    ax.plot(fitness, "k.-", alpha=0.8)
    ax.set_ylabel("Fitness (km^2/s^2)")
    ax.set_xlabel("Index of inversion")


if __name__ == "__main__":
    main()
