import os
import sys
import numpy as np
from collections import OrderedDict as odict
import matplotlib.pyplot as plt

# Internal
from CMBspectrum import CMBspectrum

try:
    import clik
except ImportError:
    raise ImportError("The likelihood code seems not to have been installed "+
                      "in your system.")

class Likelihood_Planck():
    """
    Class for calculating log-likelihoods.

    Once initialised, the method 'Likelihood_Planck.get_loglik(spectrum)' can be
    called any number of times for different 'spectrum' (instances of 'CMBspectrum').

    Mandatory arguments:
    --------------------

    base_folder: str
        Folder in which the likelihood data folders/files are found.

    Optional arguments:
    -------------------

    likelihoods: list of elements from ["commander", "camspec", "lowlike"]
        Names of the likelihoods to be computed.

    """
    def __init__(self, base_folder=None, likelihoods=None):
        fullnames = odict([["commander", "commander_v4.1_lm49.clik"],
                           ["camspec",   "CAMspec_v6.2TN_2013_02_26_dist.clik"],
                           ["lowlike",   "lowlike_v222.clik"]])
        likelihoods_fullnames = []
        if likelihoods:
            for lik in likelihoods:
                if lik.lower() in fullnames:
                    likelihoods_fullnames.append(fullnames[lik.lower()])
                else:
                    raise ValueError("Likelihood name not recognised: %s.\n"%lik+ 
                                     "Valid ones are "+str(fullnames.keys()))
            self._likelihoods_names = likelihoods_fullnames
        else:
            self._likelihoods_names = fullnames.values()
        # Initialize!
        self._likelihoods = odict()
        for lik in self._likelihoods_names:
            full_path = os.path.join(base_folder, lik)
            try:
                self._likelihoods[lik] = clik.clik(full_path)
            except clik.lkl.CError:
                raise ValueError("'clik' failed to initialise the requested "+
                                 "likelihood %s"%lik+", probably because it was"+
                                 " not found on the given folder: '%s'"%full_path)
        # Get nuisance parameters
        self._nuisance_parameters = dict([lik_name,{}]
                                         for lik_name in self._likelihoods_names)
        for lik in self._likelihoods_names:
            names = self._likelihoods[lik].extra_parameter_names
            self._nuisance_parameters[lik] = ({} if not names else
                odict([[pname,None] for pname in names]))


    # Interface methods #########################################################            
    def set_nuisance(self, n_dict=None, n_file=None):
        """
        Set the value of the nuisance parameters.
        Specify a dictionary via "n_dict" as "nuisance['param']=value"
        or a file name which contains the parameter values in different lines as
        'param = value'.
        """
        assert n_dict or n_file and not(n_dict and n_file), \
            ("A dictionary of values as 'n_dict={...}' OR a file name as"+
             +"'n_file='...' must be specified.")
        if n_file:
            nuisance_dict = {}
            try:
                nui = open(n_file, "r")
            except IOError:
                raise IOError("Nuisance parameters file not found: "+n_file)
            err_par = "Some parameter definition is not correctly formatted: "
            for line in nui:
                if line.strip() and line.strip()[0] != "#":
                    aux = [a.strip() for a in line.split()]
                    assert aux[1] == "=", (
                        "Some parameter definition is not correctly formatted:"+
                        " line: '%s'")
                    par, val = aux[0], aux[2]
                    try :
                        val = float(val)
                    except ValueError:
                        raise ValueError("Some parameter definition is not correctly formatted:"+
                                         " line: '%s'")
                    nuisance_dict[par] = val
        if n_dict:
            nuisance_dict = n_dict
        # Both cases, fill values
        for lik in self._likelihoods_names:
            for p in self._nuisance_parameters[lik]:
                try:
                    self._nuisance_parameters[lik][p] = nuisance_dict[p]
                except KeyError:
                    raise KeyError("Nuisance parameter '%s' not defined!"%p)

    def get_loglik(self, spectrum, verbose=False):
        """
        Returns a dictionary containing the contribution to the log-likelihood
        of each of the likelihoods requested.

        A summary of the information can be printed on screen using the keyword
        'verbose=True'.
        """
        spectrum_prepared = self._prepare_spectrum(spectrum)
        return self._get_loglik_internal(spectrum_prepared, verbose=verbose)

    def compare_loglik(self, test_CMBspectrum, reference_CMBspectrum,
                       delta_l=20, verbose=False):
        """
        dsaadfhadfhadfhadfhadf





        A summary of the information can be printed on screen using the keyword
        'verbose=True'.
        """
        # Prepare both spectra
        test_prepared      = self._prepare_spectrum(test_CMBspectrum)
        reference_prepared = self._prepare_spectrum(reference_CMBspectrum)
        # Go alog the multipoles and get the likelihood
        reference_loglik = \
            self._get_loglik_internal(reference_prepared, verbose=False)
        l_max = dict([lik, max(self._likelihoods[lik].lmax)]
                     for lik in self._likelihoods_names)
        n_cls = dict([name, len([int(i) for i in lik.has_cl if int(i)])]
                     for name, lik in self._likelihoods.items())
# TODO asegurarse de que esto esta mejor
        l_initials  = np.arange(0, 1+max(l_max.values()), delta_l)
        l_finals    = [l_initials[i+1]-1 for i in range(len(l_initials[:-1]))]
        l_finals   += [max(l_max.values())]
        l_intervals = zip(l_initials, l_finals)
        l_midpoints = [(ini+fin)/2. for ini, fin in l_intervals]
        if verbose:
            print ("Calculating likelihood differences along multipoles " +
                   "with Delta_l = %d"%delta_l)
            print "Progress: l =",
            sys.stdout.flush()
        loglik_differences = []
        for (l_ini, l_fin) in l_intervals:
            test_step = dict([lik, np.copy(reference_prepared[lik])]
                             for lik in reference_prepared)
            if verbose:
                print "[%d, %d] "%(l_ini, l_fin),
                sys.stdout.flush()
            this_loglik = {}
            for lik in self._likelihoods_names:
                # Don't calculate the likelihood more times than necessary
                if l_ini > l_max[lik]:
                    this_loglik[lik] = reference_loglik[lik]
                    continue
                for i_cl in range(n_cls[lik]):
                    test_step[lik][i_cl*(1+l_max[lik])+l_ini:
                                   i_cl*(1+l_max[lik])+l_fin+1] = \
                        np.copy(                                                  
                        test_prepared[lik][i_cl*(1+l_max[lik])+l_ini:
                                           i_cl*(1+l_max[lik])+l_fin+1])
                this_loglik[lik] = self._get_loglik_internal(test_step,
                                                             only=[lik])[lik]
            loglik_differences.append(dict([lik, (reference_loglik[lik]
                                                  -this_loglik[lik])]
                                           for lik in reference_loglik))
        if verbose:
            print ""
        # Prepare for plotting
        total_differences = [sum(v for v in ll.values())
                             for ll in loglik_differences]
        accum_differences = [sum(total_differences[:i+1])
                             for i in range(len(total_differences))]
        plt.figure()
        plt.plot(l_midpoints, total_differences,  color="blue")
        plt.plot(l_midpoints, accum_differences, color="red")
        plt.savefig("/tmp/test.png")
        return l_intervals, total_differences, accum_differences

    
    # Internal methods ##########################################################
    def _prepare_spectrum(self, spectrum):
        """
        Given a 'CMBspectrum' instance, prepares a dictionary of the spectra
        required by each likelihood, in the correct format to be feeded directly
        to 'clik'.

        The output is to be passed to 'Likelihood_Planck._get_loglik_internal()'.
        """
        # Check that the input is correct
        assert isinstance(spectrum, CMBspectrum), \
            "The spectrum provided must be an instance of 'CMBspectrum'."
        # Check that nuisance parameters are defined (if one is, all are)
        for lik in self._likelihoods_names:
            if self._nuisance_parameters[lik]:
                assert self._nuisance_parameters[lik].values()[0], (
                    "Nuisance parameters not yet defined! Set their values using "+
                    "'Likelihoods.set_nuisance()'.")
        # Format of Clik :  TT EE BB TE TB EB ( l = 0, 1, 2, ... !!!)
        l = list(spectrum.ll())
        pre = range(int(l[0]))
        l = np.array(pre + l)
        prepared = np.zeros([len(l), 6])
        # NOTICE that this sets C_0 = C_1 = 0
        prepared[2:, 0] = spectrum.lCl("TT", units="muK", l_prefactor=False)
        prepared[2:, 1] = spectrum.lCl("EE", units="muK", l_prefactor=False)
        prepared[2:, 3] = spectrum.lCl("TE", units="muK", l_prefactor=False)
        # Prepare the vectors for the likelihoods:
        vectors = {}
        for lik in self._likelihoods_names:
            vectors[lik] = []
            # Check enough multipoles
            l_max = self._likelihoods[lik].lmax
            assert len(l) >= max(l_max), (
                "Not enought multipoles for likelihood "+
                "'%s' : needs %d, got %d"%(lik, max(l_max), len(l)))
            # Which spectra
            which_cls = [int(i) for i in self._likelihoods[lik].has_cl]
            for i, cli in enumerate(which_cls):
                if cli:
                    vectors[lik] += prepared[:(1+l_max[i]), i].tolist()
            # Nuisance
            for par,val in self._nuisance_parameters[lik].items():
                vectors[lik].append(val)
        return vectors

    def _get_loglik_internal(self, spectrum_prepared, only=None, verbose=False):
        """
        Actually calculates the likelihood of a previously prepared spectrum,
        i.e. the output of 'Likelihood_Planck._prepare_spectrum()'.
        """
        likelihoods = self._likelihoods_names
        if only:
            assert all([lik in likelihoods for lik in only]), \
                "Likelihood not recognised: '%s'"%lik
            likelihoods = only
        loglik = {}
        for lik in likelihoods:
            if verbose:
                print "*** Computing : "+lik
            loglik[lik] = self._likelihoods[lik](spectrum_prepared[lik])
            if verbose:
                print "loglik  = ",loglik[lik]
                print "chi2eff = ",-2*loglik[lik]
        suma = sum(a[0] for a in loglik.values())
        if verbose:
            print "*** TOTAL :"
            print "loglik  = ",suma
            print "chi2eff = ",-2*suma
        return loglik

