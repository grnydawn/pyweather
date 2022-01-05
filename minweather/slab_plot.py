"""a simple slab-data plot

Usage: python slab_plot.py <slab-file>
"""

import os, sys
import matplotlib.pyplot as plt
import pyslabs

variable_names = ("dens", "umom", "wmom", "rhot")

def main():

    # open-slab file
    with pyslabs.open(sys.argv[1]) as data:

        # iterate every variables
        for name in variable_names:

            # read variables from miniweather
            for index, variable in enumerate(data.get_array(name)):

                # generate contour plot
                plt.contourf(variable)
                plt.savefig("%s_%d.png" % (name, index))


if __name__ == "__main__":
    sys.exit(main())
