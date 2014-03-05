import os
import sys

sys.path.append("../src")
from Chain import Chain
from plot_lik import plot_lik_2D

base_folder = "./chains"
chain_name = "planck_WP"
chain = Chain(os.path.join(base_folder, chain_name))

labels = {"omega_b":   r"$\omega_\mathrm{B}$",
          "omega_cdm": r"$\omega_\mathrm{C}$",
          "H0":        r"$H_0$",
          "tau_reio":  r"$\tau_\mathrm{reio}$",
          "A_s":       r"$A_s$",
          "n_s":       r"$n_s$",}
chain.set_parameter_labels(labels)
chain.plot_correlation(params = labels.keys(), save_file="correlations.png")
