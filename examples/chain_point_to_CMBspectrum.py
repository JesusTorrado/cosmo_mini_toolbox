import os
import sys

sys.path.append("../src")
from Chain import Chain
from plot_Cl_CMB import plot_Cl_CMB

base_folder = "./chains"
name = "planck_WP"
chain = Chain(os.path.join(base_folder, name))

class_folder = "/home/torradocacho/cosmo/code/class_v1.7.2"
spectrum = chain.CMBspectrum_from_point(chain.best_fit(),
                                        class_folder=class_folder,
                                        verbose=True)

plot_Cl_CMB([spectrum], title = "Best fit CMB spectrum")
