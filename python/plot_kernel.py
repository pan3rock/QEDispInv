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
import matplotlib as mpl


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "file_ker", help="filename of computed sensitivity kernel"
    )
    parser.add_argument(
        "--comp", default="vs", help="component of kernel (vp, vs, rho)"
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=int,
        default=0,
        help="mode of the dispersion curve to show",
    )
    parser.add_argument(
        "--plot_disp", action="store_true", help="show dispersion curves"
    )
    parser.add_argument(
        "--show_cb", action="store_true", help="show the colorbar"
    )
    parser.add_argument("--cmin", type=float)
    parser.add_argument("--cmax", type=float)
    parser.add_argument("--fmin", type=float)
    parser.add_argument("--fmax", type=float)
    parser.add_argument("--zmax", type=float)
    parser.add_argument("--vmax", type=float)
    parser.add_argument(
        "--unit_m", action="store_true", help="use the unit of meter"
    )
    parser.add_argument(
        "--sum", action="store_true", help="sum along the frequency"
    )
    parser.add_argument(
        "--nolast", action="store_true", help="without showing the last layer"
    )
    parser.add_argument("--savefig", help="filename of output figure")
    args = parser.parse_args()
    file_ker = args.file_ker
    comp_show = args.comp
    mode_show = args.mode
    show_sum = args.sum
    plot_dispersion = args.plot_disp
    cmin = args.cmin
    cmax = args.cmax
    fmin = args.fmin
    fmax = args.fmax
    zmax = args.zmax
    vmax = args.vmax
    unit_m = args.unit_m
    show_colorbar = args.show_cb
    nolast = args.nolast
    savefig = args.savefig

    fh5 = h5py.File(file_ker, "r")
    disp = fh5["disp"][()]
    kvp = fh5["kvp"][()]
    kvs = fh5["kvs"][()]
    krho = fh5["krho"][()]
    z = fh5["z"][()]
    fh5.close()

    if comp_show == "vs":
        var_show = kvs
        label = "kvs"
        kunit = "dimensionless"
    elif comp_show == "rho":
        var_show = krho
        label = "krho"
        kunit = r"km $\cdot$cm$^3$/(s$\cdot$ g)"
    elif comp_show == "vp":
        var_show = kvp
        label = "kvp"
        kunit = "dimensionless"
    else:
        raise ValueError("invalid comp")

    if unit_m:
        cs *= 1.0e3
        z *= 1.0e3
        unit = "m"
    else:
        unit = "km"

    if show_sum:
        plot_sum(z, var_show, disp, unit, zmax)
        return

    if nolast:
        var_show = var_show[:-1, :]
        z = z[:-1]
    else:
        var_show = np.vstack([var_show, var_show[-1, :]])
        z = np.append(z, z[-1] * 1.5)

    if vmax is None:
        vmax = np.amax(np.abs(var_show))
    vmin = -vmax

    modes = disp[:, 2].astype(int)
    idx = modes == mode_show
    disp = disp[idx]
    var_show = var_show[:, idx]
    freqs = disp[:, 0]
    cs = disp[:, 1]

    fig, ax = plt.subplots(layout="constrained")
    cmap = "seismic"
    ax.pcolormesh(
        freqs,
        z,
        var_show,
        cmap=cmap,
        shading="auto",
        vmin=vmin,
        vmax=vmax,
        alpha=0.8,
        edgecolors="none",
        antialiased=True,
        rasterized=True,
    )

    if zmax is None:
        zmax = z.max()
    if fmin is None:
        fmin = freqs.min()
    if fmax is None:
        fmax = freqs.max()

    ax.set_xlim([fmin, fmax])
    ax.set_ylim([0.0, zmax])

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel(f"Depth ({unit})")
    ax.tick_params("both")
    ax.invert_yaxis()

    ax.text(
        0.85,
        0.95,
        f"mode {mode_show:d}",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        fontsize=18,
    )

    if plot_dispersion:
        ax.tick_params("y", colors="r")
        ax.yaxis.label.set_color("r")
        ax2 = ax.twinx()
        ax2.plot(freqs, cs, "k-", alpha=0.8)
        ax2.tick_params("y", colors="k")
        ax2.set_xlim([fmin, fmax])
        if cmin is not None and cmax is not None:
            ax2.set_ylim([cmin, cmax])
        ax2.set_ylabel(f"Phase velocity ({unit}/s)")

    if savefig:
        plt.savefig(savefig, dpi=300)

    if show_colorbar:
        fig, ax = plt.subplots(figsize=(6, 0.8), layout="constrained")
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        fig.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=ax,
            orientation="horizontal",
            label=label + f" ({kunit})",
        )
        if savefig:
            fig.savefig("colorbar.{:s}".format(savefig.split(".")[-1]), dpi=300)

    plt.show()


def plot_sum(z, var_show, disp, unit, zmax):
    modes = disp[:, 2].astype(int)
    z_show = np.append(z, max(z[-1] * 1.1, 2 * z[-1] - z[-2]))

    fig, ax = plt.subplots(layout="constrained")
    modes_show = sorted(set(modes))
    for m in modes_show:
        idx = modes == m
        var = var_show[:, idx]
        var_sum = np.sum(var, axis=1)
        var_sum = np.append(var_sum, var_sum[-1])
        ax.step(var_sum, z_show, "-", alpha=0.8, label=f"mode {m}")
    ax.legend(loc="lower right")
    if zmax is None:
        ax.set_ylim([z_show.min(), z_show.max()])
    else:
        ax.set_ylim([0, zmax])
    ax.invert_yaxis()
    ax.set_ylabel(f"Depth ({unit})")
    ax.tick_params("both")
    plt.show()


if __name__ == "__main__":
    main()
