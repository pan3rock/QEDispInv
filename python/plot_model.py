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
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("file_model")
    parser.add_argument("--zmax", type=float)
    parser.add_argument("-o", "--out")
    args = parser.parse_args()
    file_model = args.file_model
    zmax = args.zmax
    file_out = args.out

    model = np.loadtxt(file_model)
    z = model[:, 1]
    fig, ax = plt.subplots(layout="constrained")
    if zmax:
        z = np.append(z, zmax)
    else:
        z = np.append(z, 2 * z[-1] - z[-2])

    ps = []
    for i in range(3):
        var = model[:, i + 2]
        var = np.append(var, var[-1])
        (p1,) = ax.step(var, z, "-", alpha=0.7, linewidth=2)
        ps.append(p1)
    ax.legend(
        ps, [r"$\rho$ (g/cm$^3$)", "Vs (km/s)", "Vp (km/s)"], loc="lower right"
    )
    if zmax:
        ax.set_ylim([0, zmax])
    else:
        ax.set_ylim([z[0], z[-1]])
    ax.set_ylabel("Depth (km)")
    ax.invert_yaxis()
    if file_out:
        fig.savefig(file_out, dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
