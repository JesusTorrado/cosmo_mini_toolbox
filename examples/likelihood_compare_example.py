import os
import sys

sys.path.append("../src")
from Chain import Chain
from Likelihood_Planck import Likelihood_Planck

# Prepare the likelihoods
likelihoods = ["commander", "CAMspec", "lowlike"]
clik_dir = "/home/torradocacho/codes/stable/planck_likelihood_1303"
lik = Likelihood_Planck(base_folder=clik_dir, likelihoods=likelihoods)

# Prepare the spectra
CHAINS = "/data/misc/torradocacho/chains"
base_folder = os.path.join(CHAINS, "historicas/aaot/gaussN/07_final")
chain = Chain(os.path.join(base_folder, "gaussN_var_c"))
# 2 best fits
best_fit_points = chain.best_fit(how_many=2)

class_folder = "/home/torradocacho/cosmo/code/class_v1.7.2_external_Pk"

# Overall best fit
best_fit_point = best_fit_points[1]
override_params = {"command": "python "+
                   os.path.join(class_folder,
                                "external_Pk/generate_Pk_from_u_gaussN.py")}
spectrum_bf1    = chain.CMBspectrum_from_point(best_fit_point,
                                               class_folder=class_folder,
                                               override_params=override_params,
                                               verbose=True)
override_params = {"P_k_ini type": "analytic_Pk",
                   "A_s": best_fit_point[chain.index_of_param("custom2", chain=True)],
                   "n_s": best_fit_point[chain.index_of_param("custom3", chain=True)]}
spectrum_bf1_0  = chain.CMBspectrum_from_point(best_fit_point,
                                               class_folder=class_folder,
                                               override_params=override_params,
                                               verbose=True)
nuisance_bf1  = chain.nuisance_file_from_point(best_fit_point, None)
lik.set_nuisance(n_dict=nuisance_bf1)
print "The original -loglik was: ",best_fit_point[5]

# Get the respective likelihoods
print "Reference: ",sum(lik.get_loglik(spectrum_bf1_0).values())
print "Test:      ",sum(lik.get_loglik(spectrum_bf1).values())

# Get the comparison
l_intervals, total_differences, accum_differences = \
    lik.compare_loglik(spectrum_bf1, spectrum_bf1_0, verbose=True, delta_l=25)

