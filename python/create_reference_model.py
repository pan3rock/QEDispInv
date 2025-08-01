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
from scipy.interpolate import make_smoothing_spline, interp1d


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("file_disp", help="filename of dispersion curves")
    parser.add_argument(
        "--vp2vs", type=float, default=1.732, help="ratio of vp and vs"
    )
    parser.add_argument(
        "-s", "--smooth", type=float, help="lamb for make_smooothing_spline"
    )
    parser.add_argument("--dmodel", help="filename of data model")
    parser.add_argument("--zmax", type=float)
    parser.add_argument(
        "-o", "--out", default="mref.txt", help="filename of output"
    )
    parser.add_argument("--savefig")
    args = parser.parse_args()
    file_disp = args.file_disp
    vp2vs = args.vp2vs
    file_datamodel = args.dmodel
    file_out = args.out
    lam = args.smooth
    zmax = args.zmax
    savefig = args.savefig

    disp = np.loadtxt(file_disp)
    modes = disp[:, 2].astype(int)
    disp = disp[modes == 0]

    f = disp[:, 0]
    c = disp[:, 1]

    wavelen = c / f
    dep = wavelen / 3.0
    ind = np.argsort(dep)
    dep = dep[ind]
    c = c[ind] * 1.1

    dep = np.insert(dep, 0, 0.0)
    vs = np.insert(c, 0, c[0])

    dep = np.append(dep, dep[-1] * 3.0 / 2.0)
    vs = np.append(vs, vs[-1])

    if lam is not None:
        spl = make_smoothing_spline(dep, vs, lam=lam)
        vs2 = spl(dep)
    else:
        vs2 = vs[:]

    vs3 = vs2

    if file_datamodel:
        dmodel = np.loadtxt(file_datamodel)
        model = create_model_with_dmodel(dep, vs3, dmodel)
        pass
    elif dep[-1] > 1.0:
        model = create_model_brocher(dep, vs3)
    else:
        model = create_model_nearsurface(dep, vs3, vp2vs)
    np.savetxt(file_out, model, fmt="%5.0f%15.8f%12.5f%12.5f%12.5f")

    fig, ax = plt.subplots(layout="constrained")
    ax.plot(vs, dep, "k.", alpha=0.8)
    ax.plot(vs3, dep, "r-", alpha=0.8, linewidth=2)
    if zmax:
        ax.set_ylim([dep[0], zmax])
    else:
        ax.set_ylim([dep[0], dep[-1]])
    ax.invert_yaxis()
    ax.set_xlabel("Vs (km/s)")
    ax.set_ylabel("Depth (km)")
    if savefig:
        fig.savefig(savefig, dpi=300)
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


def create_model_with_dmodel(dep, vs, dmodel):
    i1 = 0
    lines = []
    intp = interp1d(dep, vs)
    for i in range(1, dmodel.shape[0]):
        i2 = np.searchsorted(dep, dmodel[i, 1])
        for j in range(i1, i2):
            lines.append([0, dep[j], dmodel[i - 1, 2], vs[j], dmodel[i - 1, 4]])
        i1 = i2

        z1 = dmodel[i, 1] - 1.0e-7
        z2 = dmodel[i, 1]
        lines.append(
            [
                0,
                z1,
                dmodel[i - 1, 2],
                intp(z1),
                dmodel[i - 1, 4],
            ]
        )
        lines.append(
            [
                0,
                z2,
                dmodel[i, 2],
                intp(z2),
                dmodel[i, 4],
            ]
        )

    nl_d = dmodel.shape[0]
    for j in range(i1, dep.shape[0]):
        lines.append(
            [0, dep[j], dmodel[nl_d - 1, 2], vs[j], dmodel[nl_d - 1, 4]]
        )

    model = np.asarray(lines)
    model[:, 0] = np.arange(model.shape[0]) + 1
    return model


if __name__ == "__main__":
    main()
