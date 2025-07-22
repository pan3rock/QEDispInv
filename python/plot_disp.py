#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import argparse

params = {
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 14,
    "font.family": "serif",
}
plt.rcParams.update(params)


if __name__ == "__main__":
    msg = "plot dispersion curves"
    parser = argparse.ArgumentParser(description=msg)
    parser.add_argument(
        "file_disp", default=None, help="file of dispersion curves"
    )
    parser.add_argument("--file_ref")
    parser.add_argument(
        "--color", action="store_true", help="using different color for modes"
    )
    parser.add_argument("--ylim", nargs=2, type=float)
    parser.add_argument("--unit_m", action="store_true", help="yaxis in m")
    parser.add_argument("-o", "--out", default=None, help=" output figure name")
    parser.add_argument("--dpi", type=int, default=300, help="figure dpi")
    args = parser.parse_args()
    file_disp = args.file_disp
    file_ref = args.file_ref
    use_color = args.color
    unit_m = args.unit_m
    ylim = args.ylim
    file_out = args.out
    dpi = args.dpi

    disp = np.loadtxt(file_disp)
    modes = set(disp[:, 2].astype(int))

    if unit_m:
        km2m = 1.0e3
        unit = "m"
    else:
        km2m = 1.0
        unit = "km"

    if use_color:
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    else:
        colors = ["k"] * 5000

    fig, ax = plt.subplots(layout="constrained")
    for i, m in enumerate(modes):
        d = disp[disp[:, 2] == m]
        ax.plot(d[:, 0], d[:, 1] * km2m, "-", c=colors[i], label=str(m))

    if use_color:
        if len(modes) < 5:
            plt.legend()

    if file_ref:
        disp_ref = np.loadtxt(file_ref)
        modes = set(disp_ref[:, 2].astype(int))
        for i, m in enumerate(modes):
            d = disp_ref[disp_ref[:, 2] == m]
            ax.plot(d[:, 0], d[:, 1] * km2m, "--", c=colors[i])

    ax.set_xlim([np.min(disp[:, 0]), np.max(disp[:, 0])])
    if ylim:
        ax.set_ylim(ylim)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Phase velocity ({:s}/s)".format(unit))

    if file_out:
        fig.savefig(file_out, dpi=dpi)
    plt.show()
