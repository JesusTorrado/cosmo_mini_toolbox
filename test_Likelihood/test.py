import os

print "Before running: define inside the folder containing the likelihoods!"
exit()

clik_dir = "/home/jesus/cosmo/clik"

import sys
sys.path.append("..")
from Likelihoods import Likelihoods
from CMBspectrum import CMBspectrum

lik = Likelihoods(folder = clik_dir)

lik.set_nuisance(n_file="./nuisance.dat")

print "\n *** CLASS ******************* \n"
folder_class = "./output_CLASS"
sp_class = CMBspectrum(folder_class, code="class")
lik.get_loglik(sp_class, verbose=True)

print "\n *** CAMB ******************* \n"
folder_camb = "./output_CAMB"
sp_camb = CMBspectrum(folder_camb, prefix="test", code="camb")
lik.get_loglik(sp_camb, verbose=True)

