import os, sys

import matplotlib.pyplot as plt
import pyslabs
import numpy as np

variable_names = ("dens", "umom", "wmom", "rhot")

def main():

    with pyslabs.open(sys.argv[1]) as data:

        for name in variable_names:
            for time, variable in enumerate(data.get_array(name)):
                plt.contourf(variable)
                plt.savefig("%s_%d.png" % (name, time))


if __name__ == "__main__":
    sys.exit(main())
