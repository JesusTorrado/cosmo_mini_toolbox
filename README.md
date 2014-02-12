# `cosmo_mini_toolbox` README

Small (but growing) set of tools for manipulating cosmological data and statistical results, and for generating publication-quality plots.

Requires:
* `numpy` (minimal version tested: 1.6.1)
* `matplotlib` (minimal version tested: 1.2.1)
* For the likelihood calculation, the ESA Planck's likelihood code (see [](http://pla.esac.esa.int/pla/aio/planckProducts.html))

## Overview and examples of the modules

### Chain.py

* Document!

### plot_lik.py

Plots the [marginal, mean, profile] likelihood of the given chains with respect to the given parameters, on a grid, meaning
* Mean:     `(sum_i #_i * -loglik_i) / (sum_i #_i)`
* Marginal: `sum_i #_i`
* Profile:  `max(-loglik_i)`
being the sums over `i` extended to all chain points falling within a given cell, and being `#_i` the number of stops of the chain point `i`.

#### TODO

* Fix padding (specially negative parameters)
* Fix small issue with the best fit points

### PlanckLogLinearScale.py

An experimental implementation of the log+linear scale used to plot the CMB power spectrum by the ESA Planck team.

#### TODO:

* fix localisation in the matplotlib GUI
* Manually change the aspect ratio
* Automatic ticks
* 'change' and 'factor' parameters to be set manually, not hard coded
* Document

### plot_Cl_CMB.py

A tool to plot (absolute) comparisons between different CMB spectra.

#### Example

The code

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
                save_file = "planck_Cl_diffs.png"
               )

produces the following plot

![Different fits to the Planck CMB spectrum](planck_Cl_diffs.png)

#### TODO

* Add Planck data points

## License

**GPL3** (see the `LICENSE` file).
