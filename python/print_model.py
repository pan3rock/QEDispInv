#!/usr/bin/env python
import argparse
import numpy as np
import h5py


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("file_inv")
    parser.add_argument("--step", action="store_true")
    parser.add_argument("-o", "--out")
    args = parser.parse_args()
    file_inv = args.file_inv
    file_out = args.out
    use_step_model = args.step

    fh5 = h5py.File(file_inv, "r")
    model_mean = fh5["model_mean"][()]
    fh5.close()

    if use_step_model:
        model_mean = linear2step(model_mean)

    for i in range(model_mean.shape[0]):
        line = model_mean[i, :]
        print("{:5.0f}{:12.5f}{:12.5f}{:12.5f}{:12.5f}".format(*line))

    if file_out:
        np.savetxt(file_out, model_mean, fmt="%5.0f%12.5f%12.5f%12.5f%12.5f")


def linear2step(model_mean):
    z = model_mean[:, 1]
    nl = model_mean.shape[0]
    model2 = np.zeros([nl, 5])
    model2[0, :] = model_mean[0, :]
    for j in range(nl - 1):
        model2[j + 1, 1] = (z[j] + z[j + 1]) / 2
        for i in range(2, 5):
            model2[j + 1, i] = model_mean[j + 1, i]
    model2[:, 0] = np.arange(nl) + 1
    return model2


if __name__ == "__main__":
    main()
