###################################################
# Class for manipulating chains and getting info  #
# independently from the code that generated them #
###################################################

import os
import numpy as np

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
        codes = ["montepython", "cosmomc"]
        assert code.lower() in codes, "Code not known. Known codes are "+str(codes)
        self._code = code.lower()
        # MontePython case
        if self._code == "montepython":
            self._load_params_montepython(folder)
            self._chains = [os.path.join(folder, a) for a in os.listdir(folder)
                            if a[-4:]==".txt"] 
        # CosmoMC case
        if self._code == "cosmomc":
            self._load_params_cosmomc(folder, prefix)
            self._chains = [os.path.join(folder, a) for a in os.listdir(folder)
                            if a.startswith(prefix) and a[-4:]==".txt"] 
        # Points
        individual_chains = []
        for chain in self._chains:
            individual_chains.append(np.loadtxt(chain))
        self._points = np.concatenate(individual_chains)
        # Scaling -- MontePython
        if self._code == "montepython":
            for i, param in enumerate(self.varying_params() + self.derived_params()):
                self._points[:, i+2] *= float(self._params["parameter"][param][4])
        # Finding best fit(s) -- faster if done noe (only once!)
        maxloglik = self.get_min("mloglik")
        bf = [i for i, mloglik in enumerate(self.points("mloglik"))
              if mloglik == maxloglik]
        self._bf_points = [self.points()[i] for i in bf]

    # Load parameters
    def _load_params_montepython(self, chain_folder) :
        assert os.path.exists(chain_folder) ,\
            "The chain folder does not exist: "+chain_folder
        self._name = chain_folder.strip("/").split("/")[-1]
        logparam = open(os.path.join(chain_folder, "log.param"))
        cosmo_arguments = {}
        parameters      = {}
        path            = {}
        self._sorted_varying_params = [] # list of the (sorted) varying parameters
        self._sorted_derived_params = [] # list of the (sorted) derived parameters
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
    def _load_params_cosmomc(self, chain_folder, prefix) :
        assert os.path.exists(chain_folder) ,\
            "The chain folder does not exist: "+chain_folder
        self._name = chain_folder.strip("/").split("/")[-1]
        logparam = open(os.path.join(chain_folder, prefix+".inputparams"))
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
        with open(os.path.join(chain_folder, prefix+".paramnames"), "r") as pfile:
            for line in pfile:
                param = line.split()[0].strip()
                if param[-1] == "*":
                    self._sorted_derived_params.append(param)
                else:
                    self._sorted_varying_params.append(param)

    # Get chain data in a code independent way
    def name(self):
        return self._name
    def varying_params(self):
        return self._sorted_varying_params
    def derived_params(self):
        return self._sorted_derived_params
    def chain_files(self):
        return self._chains
    def index_of_param(self, param):
        """
        Returns the index of the given parameter in a chain point.
        """
        if param == "#":
            return 0
        elif param == "mloglik":
            return 1
        elif param in self.varying_params():
            return 2 + self.varying_params().index(param)
        elif param in self.derived_params():
            return (2 + len(self.varying_params()) +
                    self.derived_params().index(param))
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
            return self._points[:, self.index_of_param(param)]
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
        if param in self.varying_params():
            limits = [None, None]
            for i in [1, 2]:
               lim = self._params["parameter"][param][i].strip()
               if self._code == "montepython":
                   scale = float(self._params["parameter"][param][4])
                   limits[i-1] = float(lim)*scale if lim != '-1' else None
               elif self._code == "cosmomc":
                   limits[i-1] = float(lim)
            return limits
        elif param in self.derived_params():
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
            return_value = [a[self.index_of_param(param)] for a in self._bf_points]
        if not(more_than_one) and len(return_value)>1:
            print ("WARNING: Only one best fit required, but more than one present"+
                   " with the same value of the likelihood.\n"+
                   "Use the keyword 'more_than_one=True' to get them all.")
            return return_value[0]
        elif not(more_than_one) and len(return_value)==1:
            return return_value[0]
        elif more_than_one:
            return return_value
