#!/usr/bin/env python
import argparse
import h5py
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file_sfunc")
    parser.add_argument("--sign", action="store_true")
    args = parser.parse_args()
    file_sfunc = args.file_sfunc
    show_sign = args.sign

    fh5 = h5py.File(file_sfunc, "r")
    f = fh5["f"][()]
    c = fh5["c"][()]
    sfunc = fh5["sfunc"][()]
    samples = fh5["samples"][()]
    fh5.close()

    if show_sign:
        sfunc = np.sign(sfunc)

    fig, ax = plt.subplots()
    ax.axhline(0, c="k", linestyle="-", alpha=0.6)
    for i in range(samples.shape[0]):
        ax.axvline(samples[i], c="gray", alpha=0.6)
    ax.plot(c, sfunc, "r-", alpha=0.8)
    ax.set_xlim([c[0] - 0.01, c[-1] + 0.01])
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
