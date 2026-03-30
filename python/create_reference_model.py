#!/usr/bin/env python
import matplotlib.pyplot as plt
import argparse
import numpy as np
from scipy.interpolate import make_smoothing_spline, interp1d

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
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("file_disp", help="filename of dispersion curves")
    parser.add_argument("--vp2vs", type=float, default=1.732, help="ratio of vp and vs")
    parser.add_argument(
        "-s",
        "--smooth",
        type=float,
        help="Smoothing parameter (lambda) for the data fitting. "
        "0.0 means exact interpolation (passes through all points, may wiggle). "
        "Larger values (e.g., 0.01 to 1.0+) apply stronger smoothing, ignoring local noise to approach a straighter line. "
        "If not set, simple linear interpolation is used.",
    )
    parser.add_argument("--dmodel", help="filename of data model")
    parser.add_argument("--zmax", type=float)
    parser.add_argument("-o", "--out", default="mref.txt", help="filename of output")
    parser.add_argument("--savefig", help="name of the output figure file")
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
    dep_data = dep[ind]
    vs_data = c[ind] * 1.1

    # Apply global smoothing only to real observed data (Data Domain) ---
    norm_d = dep_data.max()
    norm_v = vs_data.max()
    dn = dep_data / norm_d
    vn = vs_data / norm_v

    if lam is not None:
        spl = make_smoothing_spline(dn, vn, lam=lam)
    else:
        spl = interp1d(dn, vn, kind="linear")

    # Generate smoothed intermediate model (from first to last valid depth)
    # Dense sampling with 50 points to ensure internal smoothness
    dep_mid = np.linspace(dep_data[0], dep_data[-1], 50)
    vs_mid = spl(dep_mid / norm_d) * norm_v
    vs_mid = np.maximum(vs_mid, 0.05)  # Prevent negative values

    # Explicitly construct mathematical transitions for shallow and deep parts (Transition Domain) ---
    # Shallow parabolic transition (Head) ---
    if dep_mid[0] > 0:
        z1 = dep_mid[0]
        v1 = vs_mid[0]
        # Extract smoothed slope at the first valid point
        k1 = (vs_mid[1] - vs_mid[0]) / (dep_mid[1] - dep_mid[0])
        k1 = max(k1, 0.0)

        # Limit surface velocity decrease to at most 5% below first valid point
        max_drop = 0.05 * v1
        drop = min(0.5 * k1 * z1, max_drop)
        v_surf = v1 - drop

        # Parabolic equation: v(z) = v_surf + drop * (z / z1)^2
        # Slope is 0 at z=0, ensuring smooth transition without excessively low velocities
        dep_head = np.linspace(0.0, z1, 15)[:-1]
        vs_head = v_surf + drop * (dep_head / z1) ** 2
    else:
        dep_head = np.array([])
        vs_head = np.array([])

    # Deep parabolic transition (Tail) ---
    target_zmax = zmax if zmax else dep_mid[-1] * 1.5
    if target_zmax > dep_mid[-1]:
        zN = dep_mid[-1]
        vN = vs_mid[-1]
        # Extract smoothed slope at the last valid point
        kN = (vs_mid[-1] - vs_mid[-2]) / (dep_mid[-1] - dep_mid[-2])
        kN = max(kN, 0.0)

        dz_total = target_zmax - zN
        # Limit deep velocity increase to at most 10%
        max_inc = 0.1 * vN
        tail_dv = min(0.5 * kN * dz_total, max_inc)

        # Parabolic equation: v(x) = vN + tail_dv * (2x - x^2)
        dep_tail = np.linspace(zN, target_zmax, 25)[1:]
        x_tail = (dep_tail - zN) / dz_total
        vs_tail = vN + tail_dv * (2 * x_tail - x_tail**2)
    else:
        dep_tail = np.array([])
        vs_tail = np.array([])

    dep_out = np.concatenate((dep_head, dep_mid, dep_tail))
    vs_out = np.concatenate((vs_head, vs_mid, vs_tail))

    dep = dep_out
    vs3 = vs_out

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
    ax.plot(vs_data, dep_data, "k.", alpha=0.8)
    ax.plot(vs3, dep, "r-", alpha=0.8, linewidth=2)
    if zmax:
        ax.set_ylim([dep[0], zmax])
    else:
        ax.set_ylim([dep[0], dep[-1]])
    ax.invert_yaxis()
    ax.set_xlabel("Vs (km/s)")
    ax.set_ylabel("Depth (km)")
    ax.set_title(
        "Reference Model (s={:.4f})".format(lam)
        if lam is not None
        else "Reference Model (Linear Interpolation)"
    )
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
        lines.append([0, dep[j], dmodel[nl_d - 1, 2], vs[j], dmodel[nl_d - 1, 4]])

    model = np.asarray(lines)
    model[:, 0] = np.arange(model.shape[0]) + 1
    return model


if __name__ == "__main__":
    main()
