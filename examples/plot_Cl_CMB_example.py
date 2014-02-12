import os
import sys

sys.path.append("../src")
from CMBspectrum import CMBspectrum
from plot_Cl_CMB import plot_Cl_CMB

base_folder = "./CMB_spectra"
spectrum_folders  = [os.path.join(base_folder, "planck")]
spectrum_folders += [os.path.join(base_folder, folder)
                    for folder in os.listdir(base_folder)
                    if folder.startswith("planck_")]
spectra = [CMBspectrum(folder) for folder in spectrum_folders]

plot_Cl_CMB(spectra, title = "Different fits to the Planck CMB spectrum",
#            save_file = "test.png"
            )
