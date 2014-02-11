import os
import sys

sys.path.append("..")
from Chain import Chain
from plot_lik import plot_lik_2D

folder = "./base"
chains = [Chain(folder)]
params=["H0",     "omega_b"]
labels=[r"$H_0$", r"$\omega_b$"]

plot_lik_2D("profile", chains, params=params, labels=labels, format = "-loglik",
#           transparent=False, save_file="./%s_%s.png"%(params[0], params[1])
           )
