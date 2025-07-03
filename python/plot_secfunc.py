#!/usr/bin/env python
import argparse
import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

params = {
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 14,
    "font.family": "serif",
}
plt.rcParams.update(params)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file_sfunc")
    parser.add_argument("--sign", action="store_true")
    parser.add_argument("--sample", action="store_true")
    parser.add_argument("--extra", action="store_true")
    parser.add_argument("--N", action="store_true")
    parser.add_argument("--xlim", nargs=2, type=float)
    parser.add_argument("-o", "--out", default=None, help=" output figure name")
    args = parser.parse_args()
    file_sfunc = args.file_sfunc
    show_sign = args.sign
    show_sample = args.sample
    show_extra = args.extra
    show_N = args.N
    xlim = args.xlim
    file_out = args.out

    fh5 = h5py.File(file_sfunc, "r")
    f = fh5["f"][()]
    c = fh5["c"][()]
    sfunc = fh5["sfunc"][()]
    samples = fh5["samples"][()]
    roots = fh5["roots"][()]
    samples_extra = fh5["samples_ext"][()]
    N = fh5["N"][()]
    fh5.close()

    if show_sign:
        sfunc = np.sign(sfunc)

    fig, ax = plt.subplots(layout="constrained")
    handles, labels = [], []
    ax.axhline(0, c="k", linestyle="-", alpha=0.6)
    if show_sample:
        for i in range(samples.shape[0]):
            p1 = ax.axvline(samples[i], c="k", alpha=0.6, linewidth=1)
        handles.append(p1)
        labels.append("samples")
    if show_extra:
        for i in range(samples_extra.shape[0]):
            p2 = ax.axvline(samples_extra[i], c="b", alpha=0.6, linewidth=1)
        handles.append(p2)
        labels.append("samples(extra)")
    p3 = plot_disp(ax, roots, f)
    handles.append(p3)
    labels.append("roots")

    if show_N:
        ax2 = ax.twinx()
        ax2.plot(samples, N, "-", c="tab:blue", label="N")
        ax2.set_ylabel("N", color="tab:blue")
        ax2.tick_params(axis="y", colors="tab:blue")
    ymax = np.nanmax(np.abs(sfunc))
    ax.plot(c, sfunc / ymax, "k-", alpha=0.8)
    ax.set_ylim([-1.2, 1.2])
    ax.set_yticks([-1, 0, 1])
    if xlim:
        ax.set_xlim(xlim)
    else:
        ax.set_xlim([c[0] - 0.01, c[-1] + 0.01])
    ax.legend(handles, labels, loc="lower right")
    ax.set_xlabel("Phase velocity (km/s)")
    ax.set_ylabel("Dispersion function")
    if file_out:
        fig.savefig(file_out, dpi=300)
    plt.show()


def plot_disp(ax, roots, f):
    for c in roots:
        p = ax.axvline(c, c="r", alpha=0.6)
    return p


if __name__ == "__main__":
    main()
