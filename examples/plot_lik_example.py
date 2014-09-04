import os
import sys

sys.path.append("../src")
from Chain import Chain
from plot_lik import plot_lik_2D

base_folder = "./chains"
chain = "planck_WP"
chains = [Chain(os.path.join(base_folder, chain))]
params=["H0",     "omega_b"]
labels=[r"$H_0$", r"$\omega_b$"]

modes = ["marginal", "profile"]
for mode in modes:
    plot_lik_2D(mode, chains, params=params, labels=labels, format = "-loglik",
                save="./%s_%s_%s.png"%(params[0], params[1], mode)
                )
