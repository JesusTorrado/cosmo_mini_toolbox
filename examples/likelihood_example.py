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
spectrum_name = "planck_WP"
spectrum = CMBspectrum(os.path.join(base_folder, spectrum_name))

# Prepare the nuisance parameters
base_folder = "./nuisance"
nuisance_name = "nuisance_planck_WP.dat"
lik.set_nuisance(n_file=os.path.join(base_folder, nuisance_name))

# Get the likelihood
loglik = lik.get_loglik(spectrum, verbose=True)

