#!/usr/bin/env python
import argparse
import numpy as np
import h5py


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("file_inv")
    parser.add_argument("-o", "--out")
    args = parser.parse_args()
    file_inv = args.file_inv
    file_out = args.out

    fh5 = h5py.File(file_inv, "r")
    model_mean = fh5["model_mean"][()]
    fh5.close()

    for i in range(model_mean.shape[0]):
        line = model_mean[i, :]
        print("{:5.0f}{:12.5f}{:12.5f}{:12.5f}{:12.5f}".format(*line))

    if file_out:
        np.savetxt(file_out, model_mean, fmt="%5.0f%12.5f%12.5f%12.5f%12.5f")


if __name__ == "__main__":
    main()
