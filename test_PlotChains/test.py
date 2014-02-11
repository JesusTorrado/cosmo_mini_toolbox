import os
import sys

sys.path.append("..")
from ChainTools import Chain, plot_2Dprofile

folder = "./base"
chains = [Chain(folder)]

plot_2Dprofile(chains, params=["H0", "omega_b"], labels=["$H_0$", "$\omega_b$"],
               format = "-loglik",
#               save_file = "./test.png"
               )
