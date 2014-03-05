###################################################
# Class for manipulating chains and getting info  #
# independently from the code that generated them #
###################################################

import os
import numpy as np
import matplotlib.pyplot as plt
import re

class Chain():
    """
    Class for manipulating chains and getting info from them, independently from
    the code that generated them.

    Mandatory arguments:
    --------------------

    folder: str
        Folder in which the chain is stored.
        In the case of a 'MontePython' chain, it is enough to specify the folder.
        In the case of a 'CosmoMC' chain, a 'prefix' must also be given.


    Optional arguments:
    -------------------

    prefix: str
        Prefix of the chain files of a 'CosmoMC' chain.

    code: one of ["MontePython" (default), "CosmoMC"]
        Code with which the chain was generated.

    """
    def __init__(self, folder=None, prefix=None, code="MontePython"):
        # Check input
        assert os.path.isdir(folder), \
            "The chain folder provided is not really a folder."
        self._folder = folder
        self._prefix = prefix
        codes = ["montepython", "cosmomc"]
        assert code.lower() in codes, \
            "Code not known. Known codes are " + str(codes)
        self._code = code.lower()
        # MontePython case
        if self._code == "montepython":
            self._load_params_montepython()
            self._chains = [os.path.join(self._folder, a)
                            for a in os.listdir(self._folder) if a[-4:]==".txt"]
        # CosmoMC case
        if self._code == "cosmomc":
            self._load_params_cosmomc()
            self._chains = [os.path.join(self._folder, a)
                            for a in os.listdir(self._folder)
                            if re.match(self._prefix+"_[0-9]+\.txt", a)] 
        # Points
        individual_chains = []
        for chain in self._chains:
            # Handle empty files:
            if os.path.getsize(chain) == 0:
                continue
            individual_chains.append(np.loadtxt(chain))
        self._points = np.concatenate(individual_chains)
        # Scaling -- MontePython
        if self._code == "montepython":
            for i, param in enumerate(self.varying_parameters() +
                                      self.derived_parameters()):
                self._points[:,i+2] *= float(self._params["parameter"][param][4])
        # Finding best fit(s) -- faster if done now (only once!)
        maxloglik = self.get_min("mloglik")
        bf = [i for i, mloglik in enumerate(self.points("mloglik"))
              if mloglik == maxloglik]
        self._bf_points = [self.points()[i] for i in bf]

    # Load parameters
    def _load_params_montepython(self):
        assert os.path.exists(self._folder) ,\
            "The chain folder does not exist: "+self._folder
        self._name = self._folder.strip("/").split("/")[-1]
        logparam = open(os.path.join(self._folder, "log.param"))
        cosmo_arguments = {}
        parameters      = {}
        path            = {}
        self._sorted_varying_params = [] # list of the (sorted) varying params
        self._sorted_derived_params = [] # list of the (sorted) derived params
        self._param_labels = {}          # list of parameter labels for plots
        for line in logparam:
            if not line.strip():
                continue
            if line.strip()[0] == "#":
                continue
            if "data.cosmo_arguments" in line:
                if ".update" in line:
                    continue
                left, right = line.split("=")
                key = left[1+left.find("["):left.find("]")].strip("'").strip('"')
                cosmo_arguments[key] = right.split(";")[0].strip().strip("'").strip('"')
            if "data.parameters" in line:
                left, right = line.split("=")
                key = left[1+left.find("["):left.find("]")].strip("'").strip('"')
                parameters[key] = [a.strip() for a in
                    right.split(";")[0].strip().lstrip("[").rstrip("]").split(",")]
                param_type = parameters[key][5].strip('"').strip("'")
                if param_type in ["cosmo", "nuisance"]:
                    self._sorted_varying_params.append(key)
                if param_type == "derived":
                    self._sorted_derived_params.append(key)
                self._param_labels[key] = key
            if "data.path" in line:
                left, right = line.split("=")
                if not "data.path" in left:
                    continue
                key = left[1+left.find("["):left.find("]")].strip("'").strip('"')
                path[key] = right.split(";")[0].strip()
        logparam.close()
        params = {"cosmo_argument": cosmo_arguments,
                  "parameter":      parameters,
                  "path":           path}
        self._params = params
    def _load_params_cosmomc(self) :
        assert os.path.exists(self._folder) ,\
            "The chain folder does not exist: " + self._folder
        self._name = self._prefix
        logparam = open(os.path.join(self._folder, self._prefix+".inputparams"))
        cosmo_arguments = {}
        parameters      = {}
        for line in logparam :
            if line[0] == "#" :
                continue
            # Varying parameters
            if "param[" in line and "=" in line :
                left, right = line.split("=")
                paramname = left.split("[")[-1].strip()[:-1]
                parameters[paramname] = [a.strip() for a in right.strip().split()]
            # Cosmo code arguments (i.e. fixed)
            if not("param[") in line and "=" in line :
                left, right = [a.strip() for a in line.split("=")]
                paramname = left
                cosmo_arguments[paramname] = right
        logparam.close()
        params = {"cosmo_argument": cosmo_arguments,
                  "parameter":      parameters}
        self._params = params
        # Columns in the chain files
        self._sorted_varying_params = [] # list of the (sorted) varying parameters
        self._sorted_derived_params = [] # list of the (sorted) derived parameters
        self._param_labels = {}          # list of parameter labels for plots
        with open(os.path.join(self._folder, self._prefix+".paramnames"), "r") as pfile:
            for line in pfile:
                param = line.split()[0].strip()
                if param[-1] == "*":
                    self._sorted_derived_params.append(param)
                else:
                    self._sorted_varying_params.append(param)
                self._param_labels[param] = r"%s"%line.split()[1].strip()


    # Get chain data in a code independent way
    def name(self):
        return self._name
    def varying_parameters(self):
        return self._sorted_varying_params
    def derived_parameters(self):
        return self._sorted_derived_params
    def parameters(self):
        return self.varying_parameters()+self.derived_parameters()
    def parameter_label(self, param):
        return self._param_labels[param]
    def set_parameter_labels(self, labels):
        """
        The argument 'labels' must be a dictionary of parameters (str)
        and their labels (str, possibly in LaTeX notation, e.g. r"$H_0$").

        Missing parameters in the dictionary are left to their default values.
        Unknown ones are ignored.

        Returns the list of parameter whose label was changed.
        """
        params = []
        for param, label in labels.items():
            if param in self.parameters():
                self._param_labels[param] = label
                params.append(param)
        return params
    def chain_files(self):
        return self._chains
    def index_of_param(self, param, chain=False):
        """
        Returns the index of the given parameter.

        If the keyword 'chain' is set to True (default: False) gives the index
        within a chain point row.
        """
        offset = 2 if chain else 0
        if param == "#" and chain:
            return 0
        elif param == "mloglik" and chain:
            return 1
        elif param in self.varying_parameters():
            return offset + self.varying_parameters().index(param)
        elif param in self.derived_parameters():
            return (offset + len(self.varying_parameters()) +
                    self.derived_parameters().index(param))
        else:
            raise ValueError("Unrecognized parameter: '"+str(param)+"'.")
    def points(self, param=None):
        """
        Possibilities:
        
        * param == None (or not defined) : all chain points, as rows
        * param == "#" : number of steps / (non-normalizaed) prob. of the sample
        * param == "mloglik" : minus log-likelihood
        * param == <param_name> : chain points for said parameter

        """
        if not param:
            return self._points
        else:
            return self._points[:, self.index_of_param(param, chain=True)]
    def get_min(self, param):
        """
        Gets the minimum value of the given parameter that the chain has reached.
        """
        return self.points(param=param).min()
    def get_max(self, param):
        """
        Gets the maximum value of the given parameter that the chain has reached.
        """
        return self.points(param=param).max()
    def get_limits(self, param):
        """
        Gets the limits *imposed* on the search for the given parameter
        (different from 'get_min', 'get_max').

        Output: [lower, upper], with 'None' where no limit was imposed.
        """
        if param in self.varying_parameters():
            limits = [None, None]
            for i in [1, 2]:
               lim = self._params["parameter"][param][i].strip()
               if self._code == "montepython":
                   scale = float(self._params["parameter"][param][4])
                   limits[i-1] = float(lim)*scale if lim != '-1' else None
               elif self._code == "cosmomc":
                   limits[i-1] = float(lim)
            return limits
        elif param in self.derived_parameters():
            raise ValueError("The parameter '%s' has no limits: "%param +
                             "it is a derived parameter.")
        else:
            raise ValueError("The parameter '%s' is not recognised."%param)
    def best_fit(self, param=None, more_than_one=False):
        """
        Returns the best fit point(s) of the chain, as a list of at least one point.
        (in case there is only more than one best fit point).

        A parameter can be specified to get only its best fit value.
        """
        if not param:
            return_value = self._bf_points
        else:
            return_value = [a[self.index_of_param(param, chain=True)]
                            for a in self._bf_points]
        if not(more_than_one) and len(return_value)>1:
            print ("WARNING: Only one best fit required, but more than one present"+
                   " with the same value of the likelihood.\n"+
                   "Use the keyword 'more_than_one=True' to get them all.")
            return return_value[0]
        elif not(more_than_one) and len(return_value)==1:
            return return_value[0]
        elif more_than_one:
            return return_value

    # Covariance matrix 
    def _calculate_covariance_matrix(self):
        """
        You shouldn't need to call this function, though you may want to test
        different calculation methods.
        """
        # total_steps = np.sum(self.points("#"))
        # weights = self.points("#")/total_steps
        # means = [np.sum(weights*self.points(p)) for p in
        #          self.varying_params()+self.derived_params()]
        # print means

        # For now, just read it from the .covmat file
        if self._code == "montepython":
            print ("TODO: at this point, the covmat is read from the file " +
                   "generated by MontePython's analysis routine, " +
                   "instead of calculated here!.")
            try:
                covmat = np.loadtxt(os.path.join(self._folder,
                                                 self._name + ".covmat"))
            except IOError:
                raise IOError("The '.covmat' file was not found. "+
                              "Maybe because the chain has never been analysed.")
            self._covmat = covmat
        else:
            raise NotImplementedError("Not implemented for '%s'"%self._code)
    def _assert_calculated_covmat(self):
        if not hasattr(self, "_covmat"):
            self._calculate_covariance_matrix()
    def covariance(self, param1=None, param2=None):
        """
        Returns the covariance between 'param1' and 'param2',
        or the full covariance matrix if called without arguments.
        """
        self._assert_calculated_covmat()
        if param1 and param2:
            return self._covmat[self.index_of_param(param1, chain=False),
                                self.index_of_param(param2, chain=False)]
        else:
            return self._covmat
    def variance(self, param):
        """
        Returns the variance of 'param'.
        """
        return self.covariance(param, param)
    def correlation(self, param1=None, param2=None):
        """
        Returns the correlation between 'param1' and 'param2',
        or the full correlation matrix if called without arguments.
        """
        self._assert_calculated_covmat()
        if param1 and param2:
            return (self.covariance(param1, param2) /
                    np.sqrt(self.variance(param1)*self.variance(param2)))
        else:
            corrmat = np.ones(shape=self._covmat.shape)
            for param1 in self.parameters():
                for param2 in self.parameters():
                    corrmat[self.index_of_param(param1),
                            self.index_of_param(param2)] = \
                        self.correlation(param1, param2)
            return corrmat
    def plot_correlation(self, params=None, save_file=None,
                         dpi=150, transparent=False, turn_labels=False,
                         fontsize_params=16):
        """
        Plots the correlation matrix.

        Arguments:
        ----------

        params: list of parameter names (default: None)
            If specified, plots the correlation matrix of only said parameters.

        save_file: str (default: None)
            File in which to save the plot. If not specified, the plot is "shown".

        dpi: int (default: 150)
            Resolution used if the plot is saved to a file

        transparent: bool (default: False)
            Makes the frame around the plot transparent.

        turn_labels: bool (default: False)
            Turn the parameter labels in the x axis, if they don't fit.

        fontsize_params: float (default: 16)
            Font size of the parameter labels.

        """
        if params:
            correlations = np.ones(shape=(len(params), len(params)))
            for i, param1 in enumerate(params):
                for j, param2 in enumerate(params):
                    correlations[i, j] = self.correlation(param1, param2)
        else:
            params = self.parameters()
            correlations = self.correlation()
        fig  = plt.figure()
        ax   = fig.add_subplot(111)
        imsh = ax.imshow(correlations, cmap="RdYlBu_r",
                         interpolation="nearest",# origin="lower",
                         aspect=1, zorder=0)
        imsh.set_clim(-1, 1)
        indices = [(i, j) for i in range(correlations.shape[0])
                          for j in range(correlations.shape[1])]
        for i, j in indices:
            if i != j:
                ax.annotate("%.2f"%correlations[i, j], xy=(j, i), 
                            horizontalalignment='center',
                            verticalalignment='center')
        plt.xticks(range(correlations.shape[0]),
                   [self.parameter_label(p) for p in params])
        plt.yticks(range(correlations.shape[1]),
                   [self.parameter_label(p) for p in params])
        plt.tick_params(labelsize=fontsize_params)
        if turn_labels:
            fig.autofmt_xdate()
        if not save_file:
            plt.show()
            return
        else:
            fig.frameon = transparent
            plt.savefig(save_file, transparent=transparent, dpi=dpi,
                        bbox_inches='tight', pad_inches=0.1)
            plt.close()    
            return
