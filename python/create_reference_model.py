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
from scipy.interpolate import make_smoothing_spline


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("file_disp", help="filename of dispersion curves")
    parser.add_argument(
        "--vp2vs", type=float, default=2.0, help="ratio of vp and vs"
    )
    parser.add_argument(
        "-s", "--smooth", type=float, help="lamb for make_smooothing_spline"
    )
    parser.add_argument(
        "--add", type=float, default=0, help="value added to vs"
    )
    parser.add_argument(
        "-o", "--out", default="mref.txt", help="filename of output"
    )
    args = parser.parse_args()
    file_disp = args.file_disp
    vp2vs = args.vp2vs
    file_out = args.out
    lam = args.smooth
    value_add = args.add

    disp = np.loadtxt(file_disp)
    modes = disp[:, 2].astype(int)
    disp = disp[modes == 0]

    f = disp[:, 0]
    c = disp[:, 1]

    wavelen = c / f
    dep = wavelen / 3.0
    ind = np.argsort(dep)
    dep = dep[ind]
    c = c[ind]

    dep = np.insert(dep, 0, 0.0)
    vs = np.insert(c, 0, c[0])

    dep = np.append(dep, dep[-1] * 3.0 / 2.0)
    vs = np.append(vs, vs[-1])

    if lam is not None:
        spl = make_smoothing_spline(dep, vs, lam=lam)
        vs2 = spl(dep)
    else:
        vs2 = vs[:]

    if value_add:
        vs3 = vs2 + value_add
    else:
        vs3 = vs2

    if dep[-1] > 1.0:
        model = create_model_brocher(dep, vs3)
    else:
        model = create_model_nearsurface(dep, vs3, vp2vs)
    np.savetxt(file_out, model, fmt="%5.0f%12.5f%12.5f%12.5f%12.5f")

    fig, ax = plt.subplots(layout="constrained")
    ax.plot(vs, dep, "k.", alpha=0.8)
    if value_add:
        ax.plot(vs2, dep, "b-", alpha=0.8, linewidth=2)
    ax.plot(vs3, dep, "r-", alpha=0.8, linewidth=2)
    ax.set_ylim([dep[0], dep[-1]])
    ax.invert_yaxis()
    ax.set_xlabel("Vs (km/s)")
    ax.set_ylabel("Depth (km)")
    plt.show()


def create_model_brocher(dep, vs):
    model = np.zeros([vs.shape[0], 5])
    vp = 0.9409 + 2.0947 * vs - 0.8206 * vs**2 + 0.2683 * vs**3 - 0.0251 * vs**4
    rho = (
        1.6612 * vp
        - 0.4721 * vp**2
        + 0.0671 * vp**3
        - 0.0043 * vp**4
        + 0.000106 * vp**5
    )
    model[:, 2] = rho
    model[:, 4] = vp
    model[:, 0] = np.arange(len(dep)) + 1.0
    model[:, 1] = dep
    model[:, 3] = vs
    return model


def create_model_nearsurface(dep, vs, vp2vs):
    a = -0.22374079
    b = 1.32248261
    c = 1.54840433
    rho = a * vs**2 + b * vs + c
    vp = vp2vs * vs
    model = np.zeros([vs.shape[0], 5])

    model[:, 2] = rho
    model[:, 4] = vp
    model[:, 0] = np.arange(len(dep)) + 1.0
    model[:, 1] = dep
    model[:, 3] = vs
    return model


if __name__ == "__main__":
    main()
