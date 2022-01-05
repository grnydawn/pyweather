import os, sys

import matplotlib.pyplot as plt
import pyslabs
import numpy as np

def main():

    with pyslabs.open(sys.argv[1]) as data:

        dens = data.get_array("dens")

        plt.contourf(dens[0])
        plt.savefig("dens0.png")
        plt.contourf(dens[1])
        plt.savefig("dens1.png")
        plt.contourf(dens[2])
        plt.savefig("dens2.png")


if __name__ == "__main__":
    sys.exit(main())
