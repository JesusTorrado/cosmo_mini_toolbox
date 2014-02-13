import os
import sys

sys.path.append("../src")
from CMBspectrum import CMBspectrum
from Likelihood_Planck import Likelihood_Planck

# Prepare the likelihoods
likelihoods = ["commander", "CAMspec", "lowlike"]
clik_dir = "/home/torradocacho/cosmo/data/planck/likelihood"
lik = Likelihood_Planck(base_folder=clik_dir, likelihoods=likelihoods)

# Prepare the spectrum
base_folder = "./CMB_spectra"
name_reference = "planck_WP"
name_test      = "bfB"
spectrum_reference = CMBspectrum(os.path.join(base_folder, name_reference))
spectrum_test      = CMBspectrum(os.path.join(base_folder, name_test))

# Prepare the nuisance parameters
base_folder = "./nuisance"
nuisance_name = "nuisance_planck_WP.dat"
lik.set_nuisance(n_file=os.path.join(base_folder, nuisance_name))

# Get the likelihood
loglik = lik.compare_loglik(spectrum_reference, spectrum_test, verbose=True)

